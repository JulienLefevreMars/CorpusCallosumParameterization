import nibabel as nb

if __name__ =="__main__":
	data_folder = "/home/julienlefevre/ownCloud/Documents/Recherche/Data/CorpusCallosum/isthme_du_corps_calleux/"
	subject_name = "corpus_callosum_mask_26c_201"
	filename = data_folder + subject_name + ".nii.gz"
	g = nb.load(data_folder + subject)
	print(g)
