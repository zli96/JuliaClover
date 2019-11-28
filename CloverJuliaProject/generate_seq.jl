using Random
rng= MersenneTwister();
# t - no of dna samples, m - motif, l - length of each sample in dna
function dna_generator(t,m,l)
    out = Array{String}(undef,t)
    if length(m)>l
        print("Motif can't be longer than the length of DNA\n")
    else
        gc=0;
        if 0.41*l-floor(0.41*l)>0.5
            req=floor(0.41*l)+1
        else
            req=floor(0.41*l)
        end
        for x in m
            if x=='C' || x=='G'
                gc=gc+1
            end
        end
        if gc>=req
            for i = 1:t
                tem = "";len=l-length(m)
                j=rand(big.(1:len))
                while len>=0
                    if len==j
                        tem=string(tem,m)
                    else
                        tem=string(tem,rand(rng,['A','T']))
                    end
                    len=len-1
                end
                out[i] = tem
            end
        else
           for i = 1:t
                tem = "";l1=req-gc;len=l-length(m)-l1
                while l1>0
                    tem=string(tem,rand(rng,['C','G']))
                    l1=l1-1
                end
                j=rand(big.(1:len))
                while len>=0
                    if len==j
                        tem=string(tem,m)
                    else
                        tem=string(tem,rand(rng,['A','T']))
                    end
                    len=len-1
                end
                out[i] = tem
            end
        end
        return out
    end
end

function make_fasta(input,file)
    open(file, "w") do io
        for x in 1:length(input)
            write(io, "> Generated Sequence "*string(x)*"\n")
            write(io, input[x]*"\n")
        end
    end
end

make_fasta(dna_generator(30,"ACACTGTAAACTGCTACTACTGGGCACTGTGG",400),"input_seq.txt")
make_fasta(dna_generator(70,"ACTTTGGGAACTGCTAACACCACACCGGGTTTGGCACTGACCTACTGACTGTAGGTGG",4000),"input_bg.txt")
