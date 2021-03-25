Jennet's readme

To connect to cmslpc:
>> kinit username@FNAL.GOV
The capital letters are important here, for some reason

>> ssh -Y username@cmslpc-sl7.fnal.gov

Once logged on to cmslpc, go to your nobackup area:
>> cd nobackup

Clone the repository:
>> git clone https://github.com/jennetd/scurve-derivative
>> cd scurve-derivative

Set up root
>> source setup.sh

Run the code
>> root -l methodComparison.C

