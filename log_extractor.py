import os

class log_extractor():
    '''
    Extract Information from a Gaussian log_file
    '''
    def __init__(self, log_path):
        self.log_path = log_path
        with open(self.log_path, 'r') as f:
            self.lines = f.readlines()
            self.lines = [line.strip() for line in self.lines]

    def check_exceptions(self):
        '''
        check for exceptions
        :return:
        '''
        self.exception = False # False if there are no exceptions

        #check for normal termination
        if 'Normal termination of Gaussian' not in self.lines[-1]:
            self.exception = True
        #check for negative frequencies
        if self.freqs:
            if any(freq < 0 for freq in self.freqs):
                self.exception = True
                print('Waring!!! Negative frequency found')



        return
    def extract_frequencies(self):
        freqs = []
        for line in self.lines:
            if 'Frequencies --' in line:
                split = line.split()
                freqs.extend([split[-3], split[-2], split[-1]])
        freqs = [float(freq) for freq in freqs]
        self.freqs = freqs
        return freqs



if __name__ == "__main__":
    extr = log_extractor("/home/student/j_spie17/molecular_prosthetics/gaussian/Prosthethics_all_11_07/NCBAEDUKSDLRAI-VIFPVBQESA-N_conf_0.log")
    extr.extract_frequencies()

    extr.check_exceptions()
    print(extr.exception)


    print('done')