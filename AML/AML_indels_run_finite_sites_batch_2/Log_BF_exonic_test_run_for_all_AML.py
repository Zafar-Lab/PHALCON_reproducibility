import re
import os
import pandas as pd

'''
This python file will go to each AML directory and go to Statistical_test folder and run the Run_for_statistical_test python file.

It takes the aml_data_file argument and passes it to the Run_for_statistical_test python file which in turn sends the same data argument
to both the PHALCON script and the SCIPhIN script.

In the SCIPhIN script, the cluster labels are not even required since we are just subsetting one row from the dataframe where the position matches
and we are finding liklelihood for each cell if they were not mutated (the total likelihood will be a sum of log of no-mutation likelihood across 
all cells, since we are claiming this model under the assumption of a false positive mutation happening, so we find lilkleihood of model under 
no mutation vs. PHALCON inferred likelihood when the site is mutated and we show that likelihood of PHALCON based model is more by doing fraction of
PHALCON based likelihood vs FP (SCIPhIN) based likelihood (i.e. we just subtract because of log)).
But in PHALCON one, we need the total likelihood after running PHALCON, hence in all AML files, the cluster labels are the ones which were obtained
as the final cluster labels when PHALCON was run for the first time on the whole dataset.
'''

folder = '/home/priya/Documents/priya/Log_bf_finite_sites/AML_indels_run_finite_sites_batch_2/'
sub_folders = [name for name in os.listdir(folder) if os.path.isdir(os.path.join(folder, name))]

for aml_file_name in sub_folders:
    print("File working on: ",aml_file_name)
    os.chdir(folder+aml_file_name+'/Statistical_test_exonic')
    os.system("python Run_for_statistical_test.py -D {}".format(aml_file_name))

# make a combined file for all log test runs performed
frames=[]
for aml_file_name in sub_folders:
    log_bf_test_df = pd.read_csv("/home/priya/Documents/priya/Log_bf_finite_sites/AML_indels_run_finite_sites_batch_2/"+aml_file_name+"/Statistical_test_exonic/Log BF test exonic.txt")
    frames.append(log_bf_test_df)


os.chdir(folder)
log_bf_test = pd.concat(frames,axis=0)
log_bf_test.to_csv('Log_bf_test_final_exonic.tsv',sep='\t',header=True,index=False)