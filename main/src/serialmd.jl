module Serialmd

# Azzera l'array (o matrice) a, in-place
# a = Array n-dimensionale
function set_to_zero!(a)
    z = zero(eltype(a))
    @simd for i in eachindex(a)
        @inbounds a[i] = z
    end
    return nothing
end


#= Fonte: https://en.wikipedia.org/wiki/Periodic_boundary_conditions
Vedi "Pratical Implementation" =#
# Filtra la distanza fra due componenti delle particelle i e j se il sistema è periodico
# dist = Distanza tra i e j, può essere semplicemente j-i o ri-definita, Float
# box_size = dimensioni del box di simulazione (lato), Float
# restricted = true se una particella che esce dal box rientra dal lato opposto,
#              false se la particella è libera di muoversi ma interagisce con le repliche virtuali, Boolean
# Output: Float
function periodic_dcomp(dcomp, box_size, restricted)
    if restricted
        return dcomp - round(dcomp / box_size) * box_size
    else
        return dcomp - round(dcomp * (1. / box_size)) * box_size
    end
end


function force_formula(distance, mass1, mass2, int_strength)
    # Distanza minima consentita, evita distanze nulle: a questa distanza le forze si annullano
    min_distance = 0.01
    if distance > min_distance
        return int_strength * (1/(distance*distance))
    else
        return zero(distance)
    end
end


#= PRE: Ogni punto si trova all'interno del box (se box finito) garantito dalla funzione step_update! che gestisce 
le posizioni uscenti dal box ad ogni iterazione. =#
# Find force on each particle
# 1/r^2 interactions: Is very simple but user can replace with anything they want
function find_forces!(forces, pos, vel, acc, dim, part_num, part_types, interaction_params, mass_parts, box_size, periodic, restrict)
    # Alloca variabili per i loop
    mass1 = .0f0 # Massa particella i
    mass2 = .0f0 # Massa particella k
    int_strength = .0f0 # Forza di interazione
    distance = .0f0 # Distanza fra i e k
    drag = 0.6f0 # Forza d'attrito
    force_strength = .0f0 # Forza esercitata da k su i
    
    # distance components, componenti della distanza per ciascun asse
    dcomps = Array{Float32}(undef, dim)       
    
    set_to_zero!(forces)
    
    # Per ogni particella i
    for i in 1:part_num        
        # Calcola la massa di i
        mass1 = mass_parts[part_types[i]] # La massa non veniva usata nel codice originale!
        # Per ogni altra particella k diversa da i
        @views for k in 1:part_num
            # Controlla che non sia la stessa particella oppure bounds overflow
            if (k == i || any(isinf, pos[k, :]) || any(isinf, pos[i, :])) continue end            
            mass2 = mass_parts[part_types[k]]
            # Forza d'interazione fra tipi di particelle
            int_strength = interaction_params[part_types[i], part_types[k]]           
            # Per ogni componente d calcola la distanza (in base al sistema scelto)
            for d in 1:dim
                dcomps[d] = pos[k, d] - pos[i, d]
                if periodic 
                    dcomps[d] = periodic_dcomp(dcomps[d], box_size, restrict)
                end
                distance += dcomps[d]*dcomps[d]
            end
            distance = sqrt(distance)
                       
            # Formula della forza, arbitraria                
            force_strength =  force_formula(distance, mass1, mass2, int_strength)
                        
            # Assegna le forze agenti su ciascuna componente
            forces[i, :] .+= dcomps .* force_strength
            
            # Re-inizializza a 0 per la prossima particella
            distance = .0f0 
            force_strength = .0f0
        end
        
        # Resistenza al movimento (simula attrito)
        forces[i, :] .-= 0.5 .* drag .* vel[i, :] .* abs.(vel[i, :])  # formula arbitraria        
    end
    return nothing # Modifica C++ style
end


#= Fonte: https://en.wikipedia.org/wiki/Periodic_boundary_conditions
Vedi "Pratical Implementation" =#
# Quando il sistema è periodico, vincola una coordinata a rimanere nello stesso intervallo
# p = Coordinata di una particella su un certo asse, Float
# box_size = dimensioni del box di simulazione (lato), Float
# Output: Float
function restrict_pos(p, box_size)
    return p - floor(p / box_size) * box_size
end


# Update position, velocity, and acceleration using Velocity Verlet Algorithm
# Can deal with infinite and finite systems
# For finite system, can be periodic or can reflect off walls
function step_update!(forces, pos, vel, acc, dim, part_num, part_types, mass_parts, dt, box_size, periodic, restrict)
    # Alloca variabili di loop
    mass = .0f0    
    new_acc = .0f0 
    
    # Per ogni particella
    for i = 1:part_num
        # Controllo bounds overflow, non aggiorna particelle fuori dal valore massimo        
        if (any(isinf, pos[i,:])) continue end
        # Massa della particella corrente
        mass = mass_parts[part_types[i]]
        # Per ogni componente
        for d = 1:dim                     
            # Algoritmo Velocity Verlet            
            pos[i, d] = pos[i, d] + vel[i, d]*dt + acc[i, d]*(dt*dt*0.5) # x(t+Δt) = x(t) + v(t)Δt + 1/2*a(t)(Δt)^2            
            new_acc = forces[i, d] / mass                      
            vel[i, d] = vel[i, d] + (acc[i, d] + new_acc)*(dt*0.5) # v(t+Δt) = v(t) + 1/2*(a(t)+a(t+Δt))Δt
            acc[i, d] = new_acc # a = F/m             
            
            # Gestione dei sistemi periodici
            bounce = 0.8 # 1 per rimbalzi perfettamente elastici
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
    return nothing
end


#= Salva i valori dell'istante t nella matrice trace.
IN: 
now = Matrice con i valori dell'istante corrente t
trace = Matrice con i valori per ogni istante
t = istante corrente
OUT:
nothing
====================================================#
function save!(now, trace, t)
    for r in 1:size(now, 1)
        for c in 1:size(now, 2)
            trace[t, r, c] = now[r, c]
        end
    end
    return nothing
end


#= Funzione principale di simulazione
INPUT:
nsteps = numero di iterazioni da eseguire
sinterval = intervallo di salvataggio dell'iterazione corrente
dt = durata del singolo istante di tempo
pos = matrice delle posizioni
vel = matrice delle velocità 
acc = matrice delle accelerazioni
masses = matrice delle masse
interactions = matrice delle forze di interazione
ptypes = matrice dei tipi di particella
bsize = dimensioni del box di simulazione
periodic = true se sistema periodico
restrict = true se sistema vincolato

OUTPUT:
saved_positions 
=====================================#
function dynamics_sim!(nsteps, sinterval, dt, pos, vel, acc, masses, interactions, ptypes, box_size, periodic, restrict)
    
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
        if (t-1) % sinterval == 0
            save_count += 1
            save!(pos, saved_positions, save_count)            
        end        
    end
    # Restituisci le posizioni delle particelle in ogni istante di tempo
    return saved_positions    
end


end # Fine modulo Serialmd