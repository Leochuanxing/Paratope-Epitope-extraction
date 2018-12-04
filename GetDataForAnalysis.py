import os
os.chdir("/home/leo/Documents/Database/Pipeline/Ready_2_2_0_1")
os.listdir()

import json
with open('ready_2_2_0_1', 'r') as f:
    data = json.load(f)
type(data)
len(data)
data[:6]
#with open('parepi_aa_1_2', 'r') as f:
#    parepi_aa_1_2 = json.load(f)
#with open('parepi_aa_2_1', 'r') as f:
#    parepi_aa_2_1 = json.load(f)
#with open('parepi_aa_1_1', 'r') as f:
#    parepi_aa_1_1 = json.load(f)
#with open('parepi_aa_1_1', 'r') as f:
#    parepi_aa_1_1 = json.load(f)
#with open('ac_contact', 'r') as f:
#    ac_contact = json.load(f)
#with open('sequence', 'r') as f:
#    sequence = json.load(f)
#parepi_aa_2_1
import json
with open('Homo', 'r') as f:
    Homo = json.load(f)
with open('Mouse', 'r') as f:
    Mouse = json.load(f)
import copy
data = copy.deepcopy(Homo)
data.extend(Mouse )
len(data)
data[:10]
###########################################################################



################################################################################
class Ready_for_DataAnalysis(object):
    def __init__(self, data_from_FrameConstraint_parepi_aa, contact_limit = 6, 
                 ratio = 0, testing_percentage = 0.1, saving_directory = ''):
        self.data = data_from_FrameConstraint_parepi_aa
        self.contact = contact_limit
        self.ratio = ratio
        self.percent = testing_percentage
        if saving_directory == '':
            print('You must designate a saving directory')
        else:
            self.directory = saving_directory
        
    def Select(self):
        data = self.data
        contact = self.contact
        ratio = self.ratio
        selected = []
        for parepi in data:
            if parepi[2] >= contact and parepi[2]/parepi[3] >= ratio:
                selected.append(parepi)
        self.selected = selected
        
    def Save(self):
        
        selected = self.selected
        percent = self.percent
        directory = self.directory
        
        import random
        from math import floor
        testing_number = floor(len(selected) * percent)
        hold_out_data = random.sample(selected, testing_number)

        
        ## delete the hold_out_data from the data
        for i in hold_out_data:
            selected.remove(i)
            
        
        import json
        import os
        
        os.chdir(directory)
        with open('training_data_for_DataAnalysis', 'w') as f:
            json.dump(selected, f)
        with open('testing_data_for_DataAnalysis', 'w') as f:
            json.dump(hold_out_data, f)
        
prepare = Ready_for_DataAnalysis(data, contact_limit=0, ratio=0, testing_percentage=0.1,
                                 saving_directory='/home/leo/Documents/Database/Pipeline/Ready_2_2_0_1')        
        
prepare.Select()       
prepare.Save()        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
