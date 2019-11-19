
using Statistics
function combine_score(scores)
    combinationScores = []
    n = length(scores)
    for i in 1:n # get combined score for the scenario where the motif appear in i sequences
        A = ones(i, n) # This is the dp matrix. the A_i_n element in the matrix will give us the final answer
        for x in 1:i
            for y in 1:n
                if(x > y)
                    A[x, y] = 0
                elseif(x == 1)
                    if(y == 1)
                        A[x, y] = (x * scores[y] * 1 + (y-x) * 0) / y
                    else
                        A[x, y] = (x * scores[y] * 1 + (y-x) * A[x,(y-1)]) / y
                    end
                else
                    A[x, y] = (x * scores[y] * A[(x-1),(y-1)] + (y-x) * A[x,(y-1)]) / y
                end
            end
        end
        push!(combinationScores, A[i, n])
    end
    return mean(combinationScores)
end
