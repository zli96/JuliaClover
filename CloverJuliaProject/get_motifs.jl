function Get_Single_Strand_Motifs(fileName)
    ssPFMs = []
    currentPFM = []
    motifNames = []
    open(fileName) do file
        title = ""
        while !eof(file)
            line = readline(file)
            if(line[1:1] == ">")
                if(title == "")
                    title = line[2:length(line)]
                end
                #record what has been read so far
                if(length(currentPFM) > 0)
                    push!(motifNames, title)
                    currentPFM = Normalize(currentPFM, 0.1)
                    push!(ssPFMs, currentPFM)
                    currentPFM = []
                end
                #update the current title after the old one is recorded
                title = line[2:length(line)]
            elseif(line[1:1] == "#")
                #do nothing since this line is comment
            else
                frequencies = split(line, r"\s")
                if(length(frequencies) == 4) #each line should only have 4 numbers seperated by spaces
                    row = []
                    for frequency in frequencies
                        d = parse(Int, frequency)
                        push!(row, d)
                    end
                    push!(currentPFM, row)
                else
                    println(line)
                    println("warning: this line has $(length(counts)) elements")
                end
                if(eof(file))
                    currentPFM = Normalize(currentPFM, 0.1)
                    push!(ssPFMs, currentPFM)
                    push!(motifNames, title)
                end
            end
        end
    end
    return ssPFMs, motifNames
end

function Normalize(matrix, psuedoCount)
    nrow = length(matrix)
    ncol = length(matrix[1])
    row = ones(ncol) * psuedoCount
    psuedoMatrix = []
    for i in 1:nrow
        push!(psuedoMatrix, row)
    end
    newMatrix = matrix + psuedoMatrix
    return Normalize(newMatrix)
end

function Normalize(matrix)
    normalizedMatrix = []
    for row in matrix
        row = row/sum(row)
        push!(normalizedMatrix, row)
    end
    return normalizedMatrix
end

function Get_Reverse_Complement(motif)
    reverseComplement = []
    for row in motif
        newRow = reverse(row)
        push!(reverseComplement, newRow)
    end
    reverse!(reverseComplement)
    return reverseComplement
end


function Get_Double_Strand_Motifs(singleStrandMotifs, realDoubleStrand)
    doubleStrandMotifs = []
    for motif in singleStrandMotifs
        strands = []
        temp_motif = deepcopy(motif)
        # add a column of zeros
        for row in motif
            push!(row, 0)
        end
        push!(strands, motif)
        if(realDoubleStrand)
            reverseComplement = Get_Reverse_Complement(temp_motif)
            for row in reverseComplement
                push!(row, 0)
            end
            push!(strands, reverseComplement)
        end

        push!(doubleStrandMotifs, strands)

    end
    return doubleStrandMotifs
end


#singleStrandMotifs = Get_Single_Strand_Motifs("test.txt")
#display(singleStrandMotifs)
#doubleStrandMotifs = Get_Double_Strand_Motifs(singleStrandMotifs, true)
#display(doubleStrandMotifs)
