Contact: yiming.7.chen@gmail.com or lbeutler@gmail.com

This code is for analyzing data from the first version of photometry rig (Gunaydin et al. 2014) and most of it was written in 2014. Since most of its functions are very basic, it runs the same in the newer version of matlabs (e.g tested on 2019a).The code is also for two channel rig settings (Beutler et al. 2017). 

The major role of this code is to take a specifically formatted text file that contains the organization of experimental data and turn these raw data files into a preprocessed and organized data matrix. Therefore, the exact format of txt file used as input is essential for using this code. Below is a template and example for how to construct this text file: 

[experiment name]						
'TREATMENT'	'experiment nick name'	technical_repeats				
[repeat name]		[datafilename]	[stimulation time (sec)]	[rigs sampling rate (Hz)] [laserpower]	 [prestim time window]	[poststim time windw] 	[data channel]	RETURN
[repeat name]		[datafilename]	[stimulation time (sec)]	[rigs sampling rate (Hz)] [laserpower]	 [prestim time window]	[poststim time windw] 	[data channel]	RETURN
[repeat name]		[datafilename]	[stimulation time (sec)]	[rigs sampling rate (Hz)] [laserpower]	 [prestim time window]	[poststim time windw] 	[data channel]	RETURN
...

Example:
AgRP_fasted_LeptinIP							
TREATMENT	fslt	1111111						
15050402_1027_fslt	15050402_15050403_1027_fslt	1370	250	10	300	1800	3	RETURN
15050403_1027_fslt	15050402_15050403_1027_fslt	1400	250	10	300	1800	5	RETURN
15072301_1027_fslt	15072301_15072302_1027_fslt	1371	250	10	300	1800	3	RETURN
15072302_1027_fslt	15072301_15072302_1027_fslt	1400	250	10	300	1800	5	RETURN
m14082601_r0929_fslt	m14082601_r0929_fslt	960	250	55	300	1800	5	RETURN
m14082901_r1001_fslt	m14082901_r1001_fslt	784	250	15	300	1800	5	RETURN
m14090102_r1004_fslt	m14090102_r1004_fslt	1155	250	10	300	1800	5	RETURN
