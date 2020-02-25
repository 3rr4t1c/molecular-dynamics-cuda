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