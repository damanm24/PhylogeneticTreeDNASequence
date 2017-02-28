using Requests

ncbiSearchURL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
#db="nuccore", term="user input term"
output_file = open("output.fasta", "a")

ncbiLinkURL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"

blastURL = "https://www.ncbi.nlm.nih.gov/blast/Blast.cgi"

ncbiFetchURL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
#db="nuccore", id="?", rettype="fasta", retmode="text"

#=
Search for users terms
Using esearch gather a list of ids
Using efetch using the ids check to see if its RNA data using type text
Search for related virus sequences
Download the sequences using efetch again with FASTA format
=#



function parseSearchResult(search_result)
#This function pulls out the ID's which are obtained from the genome databse
  pmid_set = Set()
  for result_line in split(search_result, "\n")
    pmid = match(r"<Id>([0-9]+)<\/Id>", result_line)
  	if pmid != nothing
  		push!(pmid_set, pmid[1])
  	end
  end
  id_string = createPMIDString(pmid_set)
  return id_string
end

function createPMIDString(pmid_set)
#From the IDs that are obtainend from parseSearchResult this list will link the IDs from genome to nuccore and add the first three results to the nuccore set
  nuccore_set = Set()
  additions = 0
  for pmid in pmid_set
  	elinkResult = readall(get(ncbiLinkURL; query = Dict("dbfrom" => "genome", "db" => "nuccore", "id" => "$pmid")))
    for result_line in split(elinkResult, "\n")
      pmid = match(r"<Id>([0-9]+)<\/Id>", result_line)
    	if pmid != nothing && additions <= 3
    		push!(nuccore_set, pmid[1])
        additions += 1
    	end
    end
  end
  #This loop generates one string with all the nuccore IDs
  id_string = ""
  for id in nuccore_set
  	id_string = "$id,$id_string"
  end
  id_string = chop(id_string)
  return id_string
end

function getSearchResult(organisms)
  #This function takes the search_terms the user provided and then searches for the organism on the database genome
  id_strings = String[]
  for search_terms in organisms
    pmids = Set()
    search_result = readall(get(ncbiSearchURL; query = Dict("db" => "genome", "term" => "$search_terms", "retmax" => 20)))
    id_string = parseSearchResult(search_result)
    push!(id_strings, id_string)
  end
  return id_strings
end

function printResults(locus, definition)
  #This function prints out the display to select the genomes to be used for the blast serach
  for (l,d) in zip(locus, definition)
    println("Locus is: $l and the definition is: $d")
  end
end

function parseFetchResult(fetch_result)
  #This function will parse the fetch_result by taking out the locus which contains the accession id and the name of the entry in the nuccore database
  locus = String[]
  definition = String[]
  for fetch_line in split(fetch_result, "\n")
    foundLocus = match(r"<GBSeq_locus>([\s\S]*?)<\/GBSeq_locus>", fetch_line)
    foundDefinition = match(r"<GBSeq_definition>([\s\S]*?)<\/GBSeq_definition>", fetch_line)
    if foundLocus != nothing
      push!(locus, foundLocus[1])
    end
    if foundDefinition != nothing
      push!(definition, foundDefinition[1])
    end
  end
  printResults(locus, definition)
end

function parseFASTAResult(fastaString)
  #this function writes the fasta string obtained to an output file
  write(output_file, fastaString)
end

function fetchResult(id_strings, retmode)
  #This function gets either the sequence data or the metadata via xml from the entry
  fetch_result = ""
  for id_string in id_strings
    if retmode == "xml"
      fetch_result = readall(get(ncbiFetchURL; query = Dict("db" => "nuccore", "id" => id_string, "retmode" => "$retmode")))
      parseFetchResult(fetch_result)
    end
    if retmode == "FASTA"
      fetch_result = readall(get(ncbiFetchURL; query = Dict("db" => "nuccore", "id" => id_string, "rettype" => "$retmode")))
      parseFASTAResult(fetch_result)
    end
  end
end

function blastSearch(accessionString)
  #This function runs a blast search with an accessionID obained from the user
  ###
  ###  DO BLASTN QUERY FOR A GIVEN QUERY_ID (ACCESSION NUMBER)
  ###
  # send query to BLASTN (technically using MegaBLAST)
  blast_put_url = "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi"
  blast_put_result = readall(post(blast_put_url; data = Dict("CMD" => "Put", "PROGRAM" => "blastn", "BLAST_PROGRAMS" => "megaBlast", "MEGABLAST" => "on", "DATABASE" => "nr", "QUERY" => "$accessionString")))


  # retrieve the job id and expected time to job completion
  rid = ""
  rtoe = ""
  for blast_put_line in split(blast_put_result, "\n")

  				rid_array  = match(r"    RID = (.*)$", blast_put_line)
  				if rid_array != nothing
  								rid = rid_array[1]
  				end

  				rtoe_array = match(r"    RTOE = (.*)$", blast_put_line)
  				if rtoe_array != nothing
  								rtoe = parse(Int64,rtoe_array[1])
  				end
  end

  println("RID ==> $rid")
  println("RTOE ==> $rtoe")

  # pause program for amount of time expected for job to take
  println("results should be ready in approximately $rtoe seconds...")
  sleep(rtoe)


  # check for status of job, if not ready, wait 5 more seconds and try again
  status = "waiting"

  while status == "waiting"

  				blast_wait_url = "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi"
  				blast_wait_result = readall(get(blast_wait_url; query = Dict("CMD" => "Get", "FORMAT_OBJECT" => "SearchInfo", "RID" => "$rid")))

  				for blast_wait_line in split(blast_wait_result, "\n")

  								# if status is waiting, pause for 5 seconds before trying again
  								if ismatch(r"\s+Status=WAITING", blast_wait_line)
  												println("...still waiting...")
  												sleep(5)
  												status = "waiting"

  								# if status is failed, exit
  								elseif ismatch(r"\s+Status=FAILED", blast_wait_line)
  												status = "failed"
  												println("something went very wrong; search failed.", blast_wait_line)
  												exit()

  								# if status is unknown, exit
  								elseif ismatch(r"\s+Status=UNKNOWN", blast_wait_line)
  												status = "failed"
  												println("something went kinda wrong; search expired.")

  								# if status is ready, change status flag
  								elseif ismatch(r"\s+Status=READY", blast_wait_line)
  												status = "ready"
  								end


  								# if status flag is set to ready, see if any results were returned
  								if status == "ready"

  												# if results found, then prepare for retrieval, otherwise exit
  												if ismatch(r"\s+ThereAreHits=yes", blast_wait_line)
  																status = "hits found"
  												elseif  ismatch(r"\s+ThereAreHits=no", blast_wait_line)
  																status = "no hits"
  																println("no hits found for search")
  																exit()
  												end
  								end

  				end

  end


  # retrieve the results from BLAST search
  println("retrieving results!")

  # only pull the hit table (no alignments), and retrieve results in plain text
  blast_result_url = "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi"
  blast_result = readall(get(blast_result_url; query = Dict("CMD" => "Get", "FORMAT_TYPE" => "Text", "ALIGNMENTS" => "0", "RID" => "$rid")))

  println(blast_result)
  accessionIDS = String[]

  for blast_result_line = split(blast_result, "\n")

  				# for each hit line, retrieve the accession id and the evalue
  				hit_line = match(r"[a-z]\|([^\|]+)\|.*\s{2}([0-9e.-]+)", blast_result_line)
  				if hit_line != nothing
  								acc_id = hit_line[1]
  								evalue = hit_line[2]

  								# only report those that have an evalue of 0.0
  								if evalue == "0.0"
  												println("related sequence ==> $acc_id (evalue = $evalue)")
                          push!(accessionIDS, acc_id)
  								end
  				end
  end
  return accessionIDS
end
#The program needs at least two organisms to generate a "meaningful tree"
organisms = String[]
print("Please give the scientific name of a virus\n")
o1 = chomp(readline(STDIN))
print("Please give another scientific name of a virus\n")
o2 = chomp(readline(STDIN))
push!(organisms, o1)
push!(organisms, o2)

o3 = ""
#Gets a complete list of organisms
while o3 != "qqq"
  print("Give me as many viruses as you want (to stop enter qqq)\n")
  o3 = chomp(readline(STDIN))
  if o3 != "qqq"
    push!(organisms, o3)
  end
end

id_strings = String[]
#ID strings will be fetched to provide the the lsit of sequences the user can choose from
id_strings = getSearchResult(organisms)
#Uses that string to grab the data in xml format
fetchResult(id_strings, "xml")


loci = String[]
specificLocus = ""
#This while loop will get an array for all the sequences that you might want to use for the blast search
while specificLocus != "qqq" || length(loci) < 2
  print("Enter the ID you want for us to use\n")
  specificLocus = chomp(readline(STDIN))
  if specificLocus != "qqq"
    push!(loci, specificLocus)
  end
end

accessionIDS = String[]

for id in loci
  accessionIDS = blastSearch(id)
  #Blast searcch returns an array of accession ids of sequences with eval of 0
  fetchResult(accessionIDS, "FASTA")
  #Those accession ids are used to download the sequence data in FASTA format
end


println("File Created and Appended")
close(output_file)
