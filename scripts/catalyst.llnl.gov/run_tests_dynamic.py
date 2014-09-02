import sys
import atexit
import time
import subprocess
import os.path

VERBOSE = True
USE_PDEBUG = False
DEBUG = False
if DEBUG:
	USE_PDEBUG = True
USE_DIMMAP = False
USE_DIMMAP_FOR_TUNE = True
NORUN = True

USE_CATALYST = True

if USE_DIMMAP:
	graph_dir = "/dimmap/"
else:
	if USE_CATALYST:
		graph_dir = "/l/ssd/"
	else:
		graph_dir = "/usr/localdisk/fusion/"

log_dir = "logs/"
executable_dir = "src/"
executable = "generate_graph_dynamic"
#executable = "generate_graph" #"run_bfs"

command_strings = []
test_count = 0

def log(s):
	if VERBOSE:
		print s
	with open(log_file_name, 'a') as f:
		f.write(s + "\n")

def init_test_dir():
	global log_dir
	global log_file_name
	global sbatch_file
	global executable

	time_stamp = str(time.time())

	while os.path.exists(log_dir+time_stamp):
		time_stamp = str(time.time())

	if DEBUG:
		log_dir += "debug/"
	else:
		log_dir += time_stamp + "/"
		os.makedirs(log_dir)

	log_file_name = log_dir + "run_tests.log"
	log("Test Motivation:")
	if len(sys.argv) == 2:
		log(str(sys.argv[1]))
	elif DEBUG:
		log("Debuging...")
	else:
		var = raw_input("Please test motivation: ")
		log(var)

	sbatch_file = log_dir + "batch.sh"

	if DEBUG:
		executable = executable_dir+executable
	else:
		cmd = ['cp', executable_dir+executable, log_dir+executable]
		p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
		while (p.poll() == None):
			pass
		executable = log_dir+executable

def generate_shell_file():
	block_start = "echo -e \"\\n\\n------------------------------------\"\n"
	block_end = "echo -e \"------------------------------------\\n\\n\"\n"
	i = 0

	if USE_CATALYST:
		slurm_options = " --clear-ssd "
	else:
		slurm_options = ""
	
	if USE_DIMMAP:
		slurm_options += "--di-mmap=" + str(10*256) + " "
	elif USE_DIMMAP_FOR_TUNE:
		slurm_options += "--di-mmap=" + str(10*256) + " "

	with open(sbatch_file, 'w') as f:
		f.write("#!/bin/bash\n")
		for cmd in command_strings:
			nodes = cmd[0]
			processes = cmd[1]
			cmd_str = cmd[2]

			if USE_PDEBUG:
				slurm_options += "-ppdebug -w catalyst322"
			elif (nodes >= 128):
				slurm_options += "-pdit128"
			elif (nodes >= 64):
				slurm_options += "-pdit64_1"
			elif (nodes >= 32):
				slurm_options += "-pdit36"


			if DEBUG:
				cmd_log_fname = log_dir+"test_%j.out"
				cmd_error_log_fname = log_dir+"test_%j.out"
			else:
				cmd_log_fname = log_dir+"test_"+str(i)+".out"
				cmd_error_log_fname = log_dir+"test_"+str(i)+".out"

			sbatch = "sbatch " + slurm_options + " -N" +str(nodes) + " -o" + cmd_log_fname + " -e" + cmd_error_log_fname + " << EOF \n"

			s = "#!/bin/sh\n"

			s += block_start + "echo Nodes: \n" + block_end
			s += "echo \"SLURM_NODELIST = \$SLURM_NODELIST \"\n"

			s += block_start + "echo Tuned Info: \n" + block_end
			s += "echo \"/proc/sys/vm/dirty_ratio = \$(cat /proc/sys/vm/dirty_ratio)\" \n"
			s += "echo \"/proc/sys/vm/dirty_background_ratio = \$(cat /proc/sys/vm/dirty_background_ratio)\" \n"
			s += "echo \"/proc/sys/vm/dirty_expire_centisecs = \$(cat /proc/sys/vm/dirty_expire_centisecs)\" \n"

			s += block_start + "echo free -m \n" + block_end
			s += "free -m \n"

			s += block_start + "echo Top 10 for memory using process \n" + block_end
			s += "ps alx  | awk '{printf (\"%d\\t%s\\n\", \\$8, \\$13)}' | sort -nr | head -10 \n"
			if USE_CATALYST:
				s += block_start + "echo df -h /l/ssd \n" + block_end
				s += "df -h -h /l/ssd  \n"
			else:
				s += block_start + "echo df -h /usr/localdisk/fusion \n" + block_end
				s += "df -h /usr/localdisk/fusion \n"

			s += block_start + "echo io-stat -m | grep md0 2>&1\n" + block_end
			s += "iostat -m | grep Device 2>&1 \n"
			s += "iostat -m | grep md0 2>&1 \n"

			s += "date \n"
			s += block_start + "echo Executable Log \n" + block_end
			s += "srun -N" +str(nodes) + " -n" + str(processes) + " " + cmd_str  + " \n"
			s += "date \n"

			s += block_start + "echo free -m \n" + block_end
			s += "free -m \n"

			if USE_CATALYST:
				s += block_start + "echo df -h /l/ssd \n" + block_end
				s += "df -h /l/ssd  \n"
			else:
				s += block_start + "echo df -h /usr/localdisk/fusion \n" + block_end
				s += "df -h /usr/localdisk/fusion \n"

			if USE_CATALYST:
				s += block_start + "echo du -sh /l/ssd/out.graph* \n" + block_end
				s += "du -sh /l/ssd/out.graph* \n"
			else:
				s += block_start + "echo du -sh /usr/localdisk/fusion/out.graph* \n" + block_end
				s += "du -sh /usr/localdisk/fusion/out.graph* \n"

			s += block_start + "echo io-stat -m | grep md0 2>&1\n" + block_end
			s += "iostat -m | grep Device 2>&1 \n"
			s += "iostat -m | grep md0 2>&1 \n"

			if USE_DIMMAP:
				s += block_start + "echo cat /proc/di-mmap-runtimeA-stats \n" + block_end
				s += "cat /proc/di-mmap-runtimeA-stats \n"

#			s += block_start + "echo dmesg \n" + block_end
#			s += "dmesg\n"

			if USE_CATALYST:
				s += block_start + "echo ls -lst /l/ssd/ \n" + block_end
				s += "ls -lst /l/ssd/\n"
			else:
				s += block_start + "echo ls -lst /usr/localdisk/fusion/ \n" + block_end
				s += "ls -lst /usr/localdisk/fusion/ \n"

			if USE_DIMMAP:
				s += block_start + "echo ls /dimmap/ \n" + block_end
				s += "ls /dimmap/\n"

			if USE_CATALYST:
				s += "rm /l/ssd/out.graph*\n"
			else:
				s += "rm /usr/localdisk/fusion/out.graph*\n"

			s += "EOF\n\n"



			f.write(sbatch + s+ "\n\n")

			i +=1



def execute_shell_file():
	if not NORUN:
		cmd = ['sh', sbatch_file]
		subprocess.call(cmd)

def add_command(nodes, processes, cmd):
	global test_count

	if not DEBUG:
		cmd_log_fname = log_dir+"test_"+str(test_count)+".out"

		log(str(test_count) + ":\t" + " ".join(cmd))

		with open(cmd_log_fname, 'w') as f:
			temp = "Test number: %d\n" %(test_count)
			f.write(temp)
			temp = "SRun Args: %s\n" %(" ".join(cmd))
			f.write(temp)
			temp = "Nodes: %d\n" %(nodes)
			f.write(temp)
			temp = "Processes: %d\n" %(processes)
			f.write(temp)
			f.write("\n")

	test_count += 1

	command_strings.append([nodes, processes, " ".join(cmd)])

def create_commands(initial_scale, scale_increments, max_scale,
	inital_nodes, node_multipler, max_nodes,
	intial_threshold, threshold_multiplier, data_type,
	low_deg_tlh_s, low_deg_tlh_e):


	graph_file = graph_dir+"out.graph"

	delete_ratio_list = [0, 5, 10, 50, 90]
	for k in delete_ratio_list:

		for i in range(low_deg_tlh_s, low_deg_tlh_e+1) :

			save_file = 0
			compare_files = 0
			test_type = "RMAT"
			chunk_size = 20
			edges_factor = 16
			scale = initial_scale
			nodes = inital_nodes
			degree_threshold = intial_threshold

			while (nodes <= max_nodes and (scale <= max_scale or max_scale == -1) ):
				processes = 1 * nodes

				cmd = [executable, test_type, str(scale), str(edges_factor), str(0), str(degree_threshold), graph_file, str(save_file), str(compare_files), str(chunk_size), data_type, str(i), str(k)]
				add_command(nodes, processes, cmd)

				nodes *= node_multipler
				scale += scale_increments
				degree_threshold *= threshold_multiplier


init_test_dir()


if DEBUG:
	#create_commands(17, 1, 17, 1, 1, 1, 1024, 1, "VC_VC")
	#create_commands(17, 1, 17, 1, 1, 1, 1024, 1, "MP_VC", 1, 10)
	#create_commands(17, 1, 17, 1, 1, 1, 1024, 1, "RB_HS")
	create_commands(14, 1, 15, 1, 1, 1, 1024, 1, "DG_AW", 1, 10)

else:
	#create_commands(17, 1, 30, 1, 1, 1, 1024, 1)
	#create_commands(24, 1, 24, 1, 1, 1, 1024, 1, "VC_VC")
	#create_commands(25, 1, 25, 1, 1, 1, 1024, 1, "MP_VC", 1, 1)
	#create_commands(22, 1, 22, 1, 1, 1, 1024, 1, "RB_HS", 1, 1)
	create_commands(24, 1, 24, 1, 1, 1, 1024, 1, "RB_MX", 1, 3)

#Data Scaling test spawning
#create_commands(29, 1, 31, 1, 1, 1, 1024, 1)


#Weak Scaling test spawning
#create_commands(20, 2, -1, 1, 2, 64, 1024, 2)

#make bash file and run it
generate_shell_file()
execute_shell_file()

log("Finished after generating %d Srun Tasks\n" %(test_count))

