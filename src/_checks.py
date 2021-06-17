import os

def _check_file_exists(f,outdir):
	"""
	Filelist is imported from blast.py return
	"""
	exist = os.path.exists(f)
	if exist == True:
		return True
	logging.debug("Blastn result file does not exist in %s" %(outdir))
	raise Exception("Blastn result file does not exist in %s" %(outdir))