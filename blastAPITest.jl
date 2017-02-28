using Requests

blastURL = "https://www.ncbi.nlm.nih.gov/blast/Blast.cgi?QUERY=AF266290.1&DATABASE=nr&PROGRAM=blastn&FORMAT_TYPE=XML&CMD=Put"
#ncbiFetchURL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

accessionString = ["AF266290.1", "AY685921.1"]

request = readall(get(blastURL))
println(request)

#seq1 = readall(get(ncbiFetchURL; query = Dict("db" => "nuccore", "id" => accessionString[1], "retmode" => "fasta", "rettype" => "text")))
#seq2 = readall(get(ncbiFetchURL; query = Dict("db" => "nuccore", "id" => accessionString[1], "retmode" => "fasta", "rettype" => "text")))
