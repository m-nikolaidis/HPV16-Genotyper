import datetime
import csv
import logging
import argparse
import sys
from Bio import SeqIO
import re
# import annot


def create_dirs(out_dir):
		"""
		Perform the necessary actions to set up the script working environment
		"""
		os.makedirs(out_dir + "/Script_out", exist_ok=False)
		os.makedirs(out_dir + "/tmp_dir", exist_ok=False)
		logging.debug(" Environment set successfully in dir %s" %(out_dir))

if __name__ == "__main__":
	now = datetime.datetime.now()
	start_dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
	print("Script started at ", start_dt_string)
	logging.info("####Starting Pipeline for Instance %s at %s #####" %(instance_name, dt_string))
	default_params = _defaultParams()
	execution_params = _updateParams(param_file, default_params)
	os.mkdir(execution_params['out'])
	logging.basicConfig(filename=execution_params['out'] + "/LOGFILE.log", format='%(asctime)s %(message)s', level=logging.DEBUG)
	# print(execution_params['in'])

	# TODO: clean tmp dirs (out_dir + "/tmpBlastn_dir/)
	# print("Cleaning tmp dirs")
	# os.rmdir()

end_dt_str = now.strftime("%d/%m/%Y %H:%M:%S")
print("The scripts started at: %s \n Finished at %s" %(start_dt_string, end_dt_str))
logging.info("The scripts started at: %s \n Finished at %s" %(start_dt_string, end_dt_str))