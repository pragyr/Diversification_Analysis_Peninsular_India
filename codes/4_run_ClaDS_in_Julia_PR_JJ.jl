using PANDA
using JLD2

# set working Directory
cd("../data/all_trees_or_super_trees")

# make a directory for storing the ClaDS results
mkdir("../ClaDS_output_Julia")

# read the filenames and store them in a vector
filenames = readdir()

# get rid of hidden files. in linux they start with a "."
filenames = filenames[findall(x -> !startswith(x,"."),filenames)]

lineage_names = SubString.(filenames,1,(length.(filenames).-4))

for i in 1:length(filenames)
	tree = load_tree(filenames[i])
	println("Currently working with the tree of - ",lineage_names[i],"\n")
	output = infer_ClaDS(tree)
	@save "../ClaDS_output_Julia/"*lineage_names[i]*"_julia_ouput" output
	save_ClaDS_in_R(output, "../ClaDS_output_Julia/"*lineage_names[i]*".Rdata")
	println("Done with ",lineage_names[i])
end
