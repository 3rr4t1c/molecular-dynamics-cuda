#= Questo modulo contiene tutte funzioni utili per gestire
un sistema di simulazione, salvataggio su disco, generazione di dati random, visualizzazione
=#
module Utilmd
using Serialization
using Random
using Plots

# Funzioni per salvare e caricare da disco tutto il necessario
# Per riavviare la simulazione

# Salva su disco
function disk_save(file_name, pos, vel, acc, masses, interactions, ptypes, box_size, periodic, restrict)
    serialize(file_name, (pos = pos,
         vel = vel, 
         acc = acc, 
         masses = masses, 
         interactions = interactions, 
         ptypes = ptypes, 
         box_size = box_size, 
         periodic = periodic, 
         restrict = restrict))
    return nothing
end

# Carica da disco
function disk_load(file_name)
    return deserialize(file_name)
end


# Generate strengths of interactions between particle types
# Very simple but user could replace with anything they wanted
# Genera le intensità delle interazioni fra tipologie di particelle
# Modifiche: Cambio tipo da Float64 a Float32, inserito break perché matrice simmetrica
function gen_interaction(num_part_types)
    interaction_params = zeros(num_part_types, num_part_types)
    rng = MersenneTwister()
    for i=1:num_part_types
        for j = 1:num_part_types
            if (i==j) # Self-interaction is randomly repulsive
                interaction_params[i, j] = -rand(rng, Float32)
            elseif (i<j) # Others randomly attractive
                val = rand(rng, Float32)
                interaction_params[i,j] = val
                interaction_params[j,i] = val
            else
                continue # Inserito continue, arresta il ciclo, matrice simmetrica
            end
        end        
    end    
    return interaction_params
end


# Genera dei dati random per avviare una simulazione
function random_data(dim, part_num, num_part_types, box_size)
    
    # Costanti dipendenti dai parametri
    interactions = gen_interaction(num_part_types) # parametri di interazione, matrice quadrata, num_part_types^2
    #masses = rand(Float32, num_part_types) .* 0.9 .+ 0.1 # Massa per ogni tipo di particella, WARNING! NO MASSE NULLE!
    masses = ceil.(rand(Float32, num_part_types) .* 3 .+ .1) # alternativa
    ptypes = rand(1:num_part_types, part_num) # Tipologia di particella per ogni particella, array di Int
      
    # Variabili aggiornate ad ogni iterazione
    pos = box_size.*rand(Float32, part_num, dim) # Initialized to be randomly placed within a box (Usare funzione?)
    vel = zeros(Float32, part_num, dim) # Initialized to zero
    acc = zeros(Float32, part_num, dim) # Initialized to zero

    return pos, vel, acc, masses, interactions, ptypes
end


### ANIMAZIONI ###

# Mappa dei colori per particella
function color_types(ptypes)
    # Lista di colori disponibili, se ci sono più tipologie di particelle ri-assegna i colori dall'inizio.
    colors = [:cyan, :lightgreen, :red, :magenta, :brown, :blue, :yellow, :black]
    colormap = Array{Symbol}(undef, size(ptypes))
    for i in eachindex(ptypes)
        colormap[i] = colors[ptypes[i] % size(colors, 1)]
    end
    return colormap
end

# Animazione per il caso tridimensionale
function gif_animate3D(saved_positions, part_types, lims, ticks, file_name)   

    pcolors = color_types(part_types)
    #msizes = round.(mass_parts3d .+ 4, digits=1)
    anim = begin 
        @animate for i in 1:size(saved_positions, 1)
            scatter(saved_positions[i,:,1], 
                    saved_positions[i,:,2], 
                    saved_positions[i,:,3], 
                    lims=lims,
                    ticks=ticks,
                    legend=false, 
                    aspect_ratio=1, 
                    markercolors=pcolors)
                    #markersize=msizes)
        end every 1
    end
    return gif(anim, "./$file_name.gif", fps=30)
end

# Animazione per il caso tridimensionale
function gif_animate2D(saved_positions, part_types, lims, ticks, file_name)   

    pcolors = color_types(part_types)
    #msizes = round.(mass_parts3d .+ 4, digits=1)
    anim = begin 
        @animate for i in 1:size(saved_positions, 1)
            scatter(saved_positions[i,:,1], 
                    saved_positions[i,:,2],
                    lims=lims,
                    ticks=ticks,
                    legend=false, 
                    aspect_ratio=1, 
                    markercolors=pcolors)
                    #markersize=msizes)
        end every 1
    end
    return gif(anim, "./$file_name.gif", fps=30)
end


end # Fine modulo Utilmd