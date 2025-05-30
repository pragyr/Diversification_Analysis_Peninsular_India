#!/usr/bin/env -S julia
using PANDA
using JLD2

# set working Directory
cd("data/trees_for_ClaDS_analyses")

# make a directory for storing the ClaDS results
if !isdir("../ClaDS_output_Julia")
	mkdir("../ClaDS_output_Julia")
end

# read the filenames and store them in a vector
filenames = readdir()

# get rid of hidden files. in linux they start with a "."
filenames = filenames[findall(x -> !startswith(x,"."),filenames)]

# remove the ".tre" from the filenames to retain the lineage/clade names
lineage_names = SubString.(filenames,1,(length.(filenames).-4))

# run the ClaDS analysis by looping over all the files
for i in 1:length(filenames)
	tree = load_tree(filenames[i]) # load the tree file
	println("Currently working with the tree of - ",lineage_names[i],"\n")
	output = infer_ClaDS(tree) # perform the ClaDS analysis
	@save "../ClaDS_output_Julia/"*lineage_names[i]*"_julia_ouput" output # save the output as a file
	save_ClaDS_in_R(output, "../ClaDS_output_Julia/"*lineage_names[i]*".Rdata") # save the output as an R project to facilitate loading the output variable 
	println("Done with ",lineage_names[i])
end
