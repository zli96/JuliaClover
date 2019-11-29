
#using Statistics
#function combine_score1(scores)
#    combinationScores = []
#    n = length(scores)
#    prev = ones(1,n)
#    current = ones(1,n)
#    for i in 1:n # get combined score for the scenario where the motif appear in i sequences
#        if(i == 1)
#            current[1] = scores[1]
#            for j in 2:n
#                current[j] = (scores[j] + (j-i)*current[j-1])/j
#            end
#        else
#            for j in 1:n
#                if(i > j)
#                    current[j] = 0
#                else
#                    current[j] = (i*scores[j]*prev[(j-1)] + (j-i)*current[j-1])/j
#                end
#            end
#        end
#        push!(combinationScores, current[n])
#        prev = deepcopy(current)
#    end
#    println("original: $combinationScores")
#    return mean(combinationScores)
#end

#a = [0.1,0.2,0.3,0.4,0.1,0.2,0.3,0.3,0.5,0.76]
#b = combine_score(a)
#println(b)
#println("!!!!!!!!!!!!!!!!!!")

function combine_score(scores)
    combinationScores = []
    n = length(scores)
    prevCol = zeros(1, n)

    for j in 1:n
        prevCell = -1
        currentCell = -1
        for i in 1:j
            if(i == 1)
                currentCell = (scores[j] + (j-i)*prevCol[i])/j
            else
                prevCell = currentCell
                currentCell = (i*scores[j]*prevCol[i-1] + (j-i)*prevCol[i])/j
                prevCol[i-1] = prevCell
            end
            #println("prevCell at[$i, $j]: $prevCell")
        end
        prevCol[j] = currentCell
    end
    return mean(prevCol)
end
