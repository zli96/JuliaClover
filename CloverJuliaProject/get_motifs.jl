include("globalVariable.jl")
function Get_Single_Strand_Motifs(fileName, pseudoCount)
    ssPFMs = []
    currentPFM = zeros(Float32, 0, ALPHSIZE)
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
                if(size(currentPFM, 1) > 0)
                    push!(motifNames, title)
                    #display(currentPFM)
                    currentPFM = Normalize(currentPFM, pseudoCount)
                    push!(ssPFMs, currentPFM)
                    currentPFM = zeros(Float32, 0, ALPHSIZE)
                end
                #update the current title after the old one is recorded
                title = line[2:length(line)]
            elseif(line[1:1] == "#")
                #do nothing since this line is comment
            else
                frequencies = split(line, r"\s")
                if(length(frequencies) == ALPHSIZE) #each line should only have ALPHSIZE numbers seperated by spaces
                    frequencies = map(x->parse(Float32, x), frequencies)
                    frequencies = reshape(frequencies, 1, ALPHSIZE)
                    currentPFM = vcat(currentPFM, frequencies)
                else
                    println(line)
                    println("warning: this line has $(length(counts)) elements")
                end
                if(eof(file))
                    currentPFM = Normalize(currentPFM, pseudoCount)
                    push!(ssPFMs, currentPFM)
                    push!(motifNames, title)
                end
            end
        end
    end
    return ssPFMs, motifNames
end

function Normalize(matrix, psuedoCount)
    nrow = size(matrix, 1)
    ncol = size(matrix, 2)
    psuedoMatrix = ones(Float32, nrow, ncol) * psuedoCount
    newMatrix = matrix + psuedoMatrix
    return Normalize(newMatrix)
end

function Normalize(matrix)
    for i in 1:size(matrix, 1)
        r = sum(matrix[i,:])
        matrix[i,:] /= r
    end
    return matrix
end

function Get_Double_Strand_Motifs(singleStrandMotifs, realDoubleStrand)
    doubleStrandMotifs = []
    for motif in singleStrandMotifs
        strands = []
        temp_motif = deepcopy(motif)
        # add a column of zeros
        nrow = size(motif, 1)
        #fifthColumn = zeros(Float64, nrow, 1)
        #motif = hcat(motif, fifthColumn)
        push!(strands, motif)
        if(realDoubleStrand)
            reverseComplement = rot180(temp_motif)
            #reverseComplement = hcat(reverseComplement, fifthColumn)
            push!(strands, reverseComplement)
        end
        push!(doubleStrandMotifs, strands)
    end
    return doubleStrandMotifs
end

#(singleStrandMotifs, motifNames) = @timeit to "ss_motif" Get_Single_Strand_Motifs("test.txt", 0.375)
##display(singleStrandMotifs[1])
#doubleStrandMotifs =  @timeit to "ds_motif" Get_Double_Strand_Motifs(singleStrandMotifs, true)
#display(doubleStrandMotifs[1][2])
