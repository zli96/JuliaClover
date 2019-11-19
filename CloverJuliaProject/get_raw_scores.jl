include("get_motifs.jl")
include("get_seqs.jl")
include("scan_seq.jl")
include("main.jl")
include("combine_score.jl")

singleStrandMotifs = Get_Single_Strand_Motifs("test.txt")
doubleStrandMotifs = Get_Double_Strand_Motifs(singleStrandMotifs, true)
#print("doubleStrandMotifs", doubleStrandMotifs)
sequenceNames = []
sequence = []
(sequence, sequenceNames) = get_seqs("fasta.txt", sequence, sequenceNames)
#print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
#print(length(sequence))
seq_info = seq_set_info(sequence)
base_probs = Vector{Float64}()
get_base_probs(sequence, base_probs)

dp=Array{Float64}(undef,length(sequence[1])+1)

#p_value_temp = Array{Float64}(undef,1)
#seq_scores_temp = Array{Float64}(undef,1)
#temp_result = result()
#print(temp_result)
results = Array{result, 1}(undef, length(doubleStrandMotifs))

motnum= 1 #temp
seqnum= 1 #temp
hit_thresh=6 #temp
base_probs=[0.3 0.1 0.5 0.1]
tot_score=scan_seq(sequence[1], doubleStrandMotifs[1], base_probs, seqnum, motnum, hit_thresh)
print("tot_score:", tot_score,"\n")

#for m in 1:length(doubleStrandMotifs)
#   for s in 1:length(sequence)
#        results[m] = result(0,Float64(0),[0.0],[0.0])
#        results[m].seq_scores = vcat(results[m].seq_scores, scan_seq(sequence[s], doubleStrandMotifs[m], base_probs, seqnum, motnum, hit_thresh))
#        results[m].raw_score = combine_scores(results[m].seq_scores)
#   end
#end
