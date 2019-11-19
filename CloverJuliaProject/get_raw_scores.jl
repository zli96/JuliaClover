include("get_motifs.jl")
include("get_seqs.jl")
include("scan_seq.jl")
include("main.jl")
include("combine_score.jl")

function print_result(sequenceFileName, seq_info)
    print("sequence file: %s (%n sequences, %n bp, %n% GC content)", sequenceFileName, seq_info.num, )
end

singleStrandMotifs = Get_Single_Strand_Motifs("test.txt")
doubleStrandMotifs = Get_Double_Strand_Motifs(singleStrandMotifs, true)
#print("doubleStrandMotifs", doubleStrandMotifs)
sequenceNames = []
sequence = []
(sequence, sequenceNames) = get_seqs("fasta.txt", sequence, sequenceNames)
#print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
#print(length(sequence))
seqInfo = init_seq_info(sequence)
print_result("fasta", seqInfo)
base_probs = Vector{Float64}()
base_probs = get_base_probs(sequence, base_probs)

dp=Array{Float64}(undef,length(sequence[1])+1)

results = Vector{result}(undef,length(doubleStrandMotifs))
# print(methods(result))
motnum= 1 #temp
seqnum= 1 #temp
hit_thresh=6 #temp

print("\nlength of ds: ", length(doubleStrandMotifs))
for m in 1:length(doubleStrandMotifs)
    sequenceScores = []
    for s in 1:length(sequence)
        push!(sequenceScores, scan_seq(sequence[s], doubleStrandMotifs[m], base_probs, seqnum, motnum, hit_thresh))
    end
    rawScore = combine_score(sequenceScores)
    results[m] = result(0, rawScore, [], sequenceScores)
end
bg_info=[]
bg_seqs=[]
junk=[]
(bg_seqs, junk) = get_seqs("bg_seqs.txt", bg_seqs, junk)
bg_info = vcat(bg_info,seq_set_info(bg_seqs))
p_vals=shuffle_bgseq(sequence,bg_seqs,doubleStrandMotifs)

# get_hits(sequence,doubleStrandMotifs,base_probs,p_vals)
