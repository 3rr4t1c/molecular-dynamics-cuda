module CPUMTmd

using Base.Threads
using Random


#= Funzione ausiliaria per verificare se una posizione non è valida.
In alcuni sistemi se la particella supera il valore maggiore possibile
per una delle sue componenti, questo viene settato ad Inf, 
andando a "contagiare" tutte le altre particelle,
risultando in una serie di valori NaN.
==============================#
function anyinfpos(i, pos, dim)

    for d in 1:dim
        if isinf(pos[i, d])
            return true
        end
    end
    return false

end


#= Fonte: https://en.wikipedia.org/wiki/Periodic_boundary_conditions, vedi "Pratical Implementation"
Filtra la distanza fra due componenti delle particelle i e j se il sistema è periodico
dist = Distanza tra i e j, può essere semplicemente j-i o ri-definita
box_size = dimensioni del box di simulazione (lato)
restricted = true se una particella che esce dal box rientra dal lato opposto
false se la particella è libera di muoversi ma interagisce con le repliche virtuali
===================================================#
function periodic_dcomp(dcomp, box_size, restricted)

    if restricted
        return @fastmath dcomp - round(dcomp / box_size) * box_size
    else
        return @fastmath dcomp - round(dcomp * (1.0f0 / box_size)) * box_size
    end

end


#= Formula della forza, può essere reimplementata a piacimento
Viene utilizzata in find_forces!
===========================================================#
function force_formula(distance, mass1, mass2, int_strength, min_distance)
    
    if distance > min_distance
        return @fastmath int_strength * (1.0f0 / (distance*distance))
    else
        return .0f0
    end

end


#= PRE: Ogni punto si trova all'interno del box, questo è garantito dalla funzione step_update! 
che gestisce le posizioni uscenti dal box ad ogni iterazione. 
Calcola la forza agente su ogni particella (eccetto quelle in posizioni non valide)
====================================================================================================================================#
function find_forces!(forces, pos, vel, acc, dim, part_num, part_types, interaction_params, mass_parts, box_size, periodic, restrict)
    
    # Azzera forze
    forces .= .0f0
    
    # Costanti
    drag = .6f0 # Coefficiente d'attrito
    min_distance = .01f0 # Distanza minima consentita

    # Per ogni particella i
    @threads for i in 1:part_num
        
        # Alloca variabili per il thread corrente
        mass1 = .0f0 # Massa particella i
        mass2 = .0f0 # Massa particella k
        int_strength = .0f0 # Forza di interazione
        distance = .0f0 # Distanza fra i e k        
        force_strength = .0f0 # Forza esercitata da k su i

        @inbounds begin
            # Calcola la massa di i
            mass1 = mass_parts[part_types[i]] # La massa non veniva usata nel codice originale!
            # Per ogni altra particella k diversa da i
            for k in 1:part_num
                # Controlla che non sia la stessa particella oppure bounds overflow
                if (k == i || anyinfpos(i, pos, dim) || anyinfpos(k, pos, dim) ) continue end            
                mass2 = mass_parts[part_types[k]]
                # Forza d'interazione fra tipi di particelle
                int_strength = interaction_params[part_types[i], part_types[k]]           
                # Per ogni componente d calcola la distanza (in base al sistema scelto)
                for d in 1:dim
                    dcomp = pos[k, d] - pos[i, d]
                    if periodic 
                        dcomp = periodic_dcomp(dcomp, box_size, restrict)
                    end
                    distance += dcomp*dcomp
                end
                distance = @fastmath sqrt(distance)
                        
                # Formula della forza, arbitraria                
                force_strength = force_formula(distance, mass1, mass2, int_strength, min_distance)
                            
                # Assegna le forze agenti su ciascuna componente                
                for d in 1:dim
                    dcomp = pos[k, d] - pos[i, d]
                    if periodic 
                        dcomp = periodic_dcomp(dcomp, box_size, restrict)
                    end
                    forces[i, d] += dcomp * force_strength
                end
                
                # Re-inizializza a 0 per la prossima particella
                distance = .0f0 

            end
            
            # Resistenza al movimento (simula attrito), formula arbitraria            
            for d in 1:dim
                forces[i, d] -= @fastmath .5f0 * drag * vel[i, d] * abs(vel[i, d])
            end

        end
    end
    return nothing

end


#= Fonte: https://en.wikipedia.org/wiki/Periodic_boundary_conditions, vedi "Pratical Implementation"
# Quando il sistema è periodico, vincola una componente a rimanere nello stesso intervallo
# p = componente di una particella su un certo asse
# box_size = dimensioni del box di simulazione (lato)
=================================#
function restrict_pos(p, box_size)

    return @fastmath p - floor(p / box_size) * box_size

end


#= Aggiorna posizione, velocità ed accelerazione usando il Velocity Verlet Algorithm.
Ogni sistema è considerato finito, per modellare un sistema infinito basta impostare restrict a false.
====================================================================================================================#
function step_update!(forces, pos, vel, acc, dim, part_num, part_types, mass_parts, dt, box_size, periodic, restrict)
    
    # Costanti
    bounce = .8f0 # 1 per rimbalzi perfettamente elastici
    
    # Alloca variabili di loop
    mass = .0f0    
    new_acc = .0f0 
    
    # Per ogni particella avvia un thread
    @threads for i in 1:part_num
        @inbounds begin
            # Controllo bounds overflow, salta le particelle che sono fuori valore massimo        
            if anyinfpos(i, pos, dim) continue end
            # Massa della particella corrente
            mass = mass_parts[part_types[i]]
            # Per ogni componente
            for d = 1:dim                     
                # Algoritmo Velocity Verlet            
                pos[i, d] = @fastmath pos[i, d] + vel[i, d]*dt + acc[i, d]*(dt*dt*.5f0) # x(t+Δt) = x(t) + v(t)Δt + 1/2*a(t)(Δt)^2            
                new_acc = forces[i, d] / mass                      
                vel[i, d] = @fastmath vel[i, d] + (acc[i, d] + new_acc)*(dt*.5f0) # v(t+Δt) = v(t) + 1/2*(a(t)+a(t+Δt))Δt
                acc[i, d] = new_acc # a = F/m             
                
                # Gestione dei sistemi periodici                
                if periodic && restrict # Tipo "Pac-Man"
                    pos[i, d] = restrict_pos(pos[i, d], box_size)                     
                elseif !periodic && restrict # Rimbalza sulle pareti
                    if pos[i, d] < 0 || pos[i, d] > box_size
                        pos[i, d] = box_size - restrict_pos(pos[i, d], box_size)                    
                        vel[i, d] = -vel[i, d] * bounce
                        acc[i, d] = -acc[i, d] * bounce
                    end                
                end           
            end
        end
    end    
    return nothing

end


#= Funzione principale di simulazione.
nsteps = numero di iterazioni da eseguire
sinterval = intervallo di salvataggio dell'iterazione corrente
track = true se si vogliono salvare le posizioni (utile per benchmark)
dt = durata del singolo istante di tempo
pos = matrice delle posizioni
vel = matrice delle velocità 
acc = matrice delle accelerazioni
masses = matrice delle masse
interactions = matrice delle forze di interazione
ptypes = matrice dei tipi di particella
box_size = dimensioni del box di simulazione
periodic = true se sistema periodico
restrict = true se sistema vincolato
RETURN:
saved_positions = posizioni in ogni istante di tempo
=====================================================#
function dynamics_sim!(nsteps::Int, 
                       sinterval::Int,
                       track::Bool, 
                       dt::Float32, 
                       pos::Array{Float32,2},
                       vel::Array{Float32,2},
                       acc::Array{Float32,2},
                       masses::Array{Float32,1},
                       interactions::Array{Float32,2},
                       ptypes::Array{Int,1},
                       box_size::Float32,
                       periodic::Bool,
                       restrict::Bool)
    
    # Determina il numero di particelle e la dimensionalità dello spazio
    part_num, dim = size(pos)
        
    # Variabili per il salvataggio
    saves_num = nsteps ÷ sinterval # Numero di posizioni da salvare
    saved_positions = Array{Float32}(undef, saves_num, part_num, dim)
    
    # Le forze da ri-calcolare in ogni istante
    forces = Array{Float32}(undef, part_num, dim) # Viene azzerato all'inizio di ogni ciclo
    
    # contatore salvataggi
    save_count = 0    
    
    # Loop per ciascun istante di tempo
    for t in 1:nsteps
        # Calcola forze agenti su ciascuna particella, (todo: parametri superflui part_num, part_types?)
        find_forces!(forces, pos, vel, acc, dim, part_num, ptypes, interactions, masses, box_size, periodic, restrict)
        # aggiornamento step
        step_update!(forces, pos, vel, acc, dim, part_num, ptypes, masses, dt, box_size, periodic, restrict)
        # Salva posizioni
        if track && (t-1) % sinterval == 0
            save_count += 1
            saved_positions[save_count,:,:] .= pos            
        end        
    end
    # Restituisci le posizioni delle particelle in ogni istante di tempo
    return saved_positions

end


#= Genera le intensità di interazione fra tipologie di particelle
num_part_types = numero di tipi di particelle
RETURN:
interaction_params = matrice dei parametri di interazione
=======================================#
function gen_interaction(num_part_types)

    interaction_params = zeros(Float32, num_part_types, num_part_types)
    rng = MersenneTwister()
    for i in 1:num_part_types
        for j in 1:num_part_types
            if (i==j) # Self-interaction is randomly repulsive
                interaction_params[i,j] = -rand(rng)
            elseif (i<j) # Others randomly attractive
                val = rand(rng)
                interaction_params[i,j] = val
                interaction_params[j,i] = val
            else
                continue # Matrice simmetrica
            end
        end        
    end    
    return interaction_params

end


#= Genera dei dati random per avviare una simulazione.
dim = numero di dimensioni
part_num = numero di particelle
num_part_types = numero di tipologie di particelle
box_size = dimensioni del box di simulazione iniziale
RETURN:
pos = matrice delle posizioni
vel = matrice delle velocità
acc = matrice delle accelerazioni
masses = vettore delle masse
interaction = matrice dei parametri di interazione
ptypes = vettore del tipo di particella 
============================================================#
function random_data(dim, part_num, num_part_types, box_size)
    
    # Costanti dipendenti dai parametri
    interactions = gen_interaction(num_part_types) # parametri di interazione, matrice quadrata, num_part_types^2
    #masses = rand(Float32, num_part_types) .* 0.9 .+ 0.1 # Massa per ogni tipo di particella, WARNING! NO MASSE NULLE!
    masses = ceil.(rand(Float32, num_part_types) .* 3 .+ .1f0) # alternativa
    ptypes = rand(1:num_part_types, part_num) # Tipologia di particella per ogni particella, array di Int
      
    # Variabili aggiornate ad ogni iterazione
    pos = box_size .* rand(Float32, part_num, dim) # Inizializza con posizioni casuali nel box 
    vel = zeros(Float32, part_num, dim) # Initialized to zero
    acc = zeros(Float32, part_num, dim) # Initialized to zero

    return pos, vel, acc, masses, interactions, ptypes

end


end # Fine modulo CPUMTmd