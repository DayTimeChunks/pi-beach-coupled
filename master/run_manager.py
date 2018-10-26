import math
import os
import shutil
import re
from copy import deepcopy
# Read in the file
lisem_runs = [183, 188, 199, 214, 225,
              242, 243, 247, 248,
              260, 269]
event_type = [2, 1, 1, 2, 1,
              3, 4, 3, 4,
              1, 1]

# Event types:
# 1 = 1 event in one day (no day break)
# 2 = event split into 2 days
# 3 = First of two consecutive days of type 1
# 4 = Second of two consecutive days of type 1
# 5 = 2 or more events in 1 day
""" Changes needed as a result of .run file changes. """
# Preparing inputs for LISEM "run" file
target_time = [162., 110., 300., 120., 180.,
               150., 270., 180., 190.,
               118., 240.]
end_time = [202., 110., 300., 300., 180.,
            150., 270., 180., 190.,
            118., 240.]

# Long time-steps only for testing!
# ts_new = [126., 126., 200., 126., 200.,
#           200., 200., 200., 200.,
#           200., 200.]

ts_new = [10., 10., 10., 10., 10.,
          10., 10., 10., 10.,
          10., 10.]
# ts_new = [9., 9., 9., 9., 9.]

assert len(lisem_runs) == \
       len(event_type) == len(end_time) ==\
       len(target_time) == len(ts_new)

interval = 99.
map_list = []
hydro_row = []

# Following loop inputs new values to "run" file AND
# tells BEACH which (i) infiltration maps and (ii) runoff value
# to get from hydro.txt row
for t in lisem_runs:
    i = lisem_runs.index(t)
    # Infiltration Maps
    elapsed_min = ts_new[i] / 60. * interval
    map_num = int(round(target_time[i] / elapsed_min))
    if target_time[i] == end_time[i]:
        map_list.append("\\infiltration.map")
    else:
        if map_num < 10:
            map_list.append("\\infilLM0.00" + str(map_num))
        else:
            map_list.append("\\infilLM0.0" + str(map_num))

    # Runoff rows
    row_ts = int(target_time[i] / (ts_new[i] / 60.))
    if target_time[i] == end_time[i]:
        hydro_row.append(-1)
    else:
        hydro_row.append(row_ts)


def create_dirs(mc_samples):
    # Clean
    old_beach_directory = os.getcwd() + "\\tss\\output"
    old_lisem_directory = os.getcwd() + "\\visuals_opL"
    if os.path.exists(old_beach_directory):
        shutil.rmtree(old_beach_directory)
    # Clean
    if os.path.exists(old_lisem_directory):
        shutil.rmtree(old_lisem_directory)

    for sample_nr in range(1, mc_samples+1):
        for ev in range(len(lisem_runs)):
            directory = os.getcwd() + "\\res\\" + str(ev)

            """ LISEM results (temporary) """
            # Clean - Delete all previous results!
            if os.path.exists(directory):
                shutil.rmtree(directory)
                shutil.copytree(os.getcwd() + "\\res\\dummy", directory)

            # For first-time-ever events, create folder repo.
            if not os.path.exists(directory):
                shutil.copytree(os.getcwd() + "\\res\\dummy", directory)

            """ BEACH MC - samples, and respective LISEM hydro.txt for each event in that sample """
            # Check tss/output/Bs folders exist.
            beach_directory = os.getcwd() + "\\tss\\output\\Bs" + str(sample_nr) + "\\" + str(ev)
            lisem_directory = os.getcwd() + "\\visuals_opL\\Ls" + str(sample_nr) + "\\" + str(ev)

            if not os.path.exists(beach_directory):
                shutil.copytree(os.getcwd() + "\\res\\dummy", beach_directory)
            if not os.path.exists(lisem_directory):
                shutil.copytree(os.getcwd() + "\\res\\dummy", lisem_directory)

            # Extra BEACH tss/output folder
            if ev == len(lisem_runs) - 1:
                beach_directory = os.getcwd() + "\\tss\\output\\Bs" + str(sample_nr) + "\\" + str(ev+1)
                # lisem_directory = os.getcwd() + "\\visuals_opL\\Ls" + str(sample_nr) + "\\" + str(ev+1)
                if not os.path.exists(beach_directory):
                    print("Created: " + beach_directory)
                    shutil.copytree(os.getcwd() + "\\res\\dummy", beach_directory)
                # if not os.path.exists(lisem_directory):
                #     print("Created: " + lisem_directory)
                #     shutil.copytree(os.getcwd() + "\\res\\dummy", lisem_directory)


def update_runfile_directories(events):
    """
    Changes the paths of all LISEM run files to the current working directory
    This is useful when making a copy of the parent module (e.g. parent: b2l_MC)

    To do so, just run this file as main (i.e. not as a module imported by beach_control).
    :return:
    """
    print(os.getcwd())

    """
    Update run files
    """
    for ev in events:
        work_dir = os.getcwd()
        new_dir = "/".join(work_dir.split("\\"))
        map_dir_end = "/maps/lisem/" + str(events.index(ev))
        rain_dir_end = "/events/"
        res_dir_end = "/res/" + str(events.index(ev))

        # Prepare LISEM initial map's folders (if they don't already exists)
        lisem_initial_dir = os.getcwd() + map_dir_end
        if not os.path.exists(lisem_initial_dir):
            shutil.copytree(os.getcwd() + "/maps/lisem/0", lisem_initial_dir)

        if event_type[events.index(ev)] == 2:
            # Prints timestep-specific map (e.g. dynamic pcraster format)
            map_zeroes = "0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0"
        else:
            # Prints only the last cumulative map
            map_zeroes = "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"

        with open('Alteck' + str(ev) + '.run', 'rt') as file:
            for old_line in file:
                work = re.match(r'Work Directory=', old_line)
                if work is not None:
                    print(old_line)
                    work_line = old_line

                map = re.match(r'Map Directory=', old_line)
                if map is not None:
                    map_line = old_line

                rain = re.match(r'Rainfall Directory=', old_line)
                if rain is not None:
                    rain_line = old_line

                rain_f = re.match(r'Rainfall file=', old_line)
                if rain_f is not None:
                    rain_f_line = old_line

                res = re.match(r'Result Directory=', old_line)
                if res is not None:
                    res_line = old_line

                end_time_txt = re.match(r'End time=', old_line)
                if end_time_txt is not None:
                    end_time_line = old_line

                time = re.match(r'Timestep=', old_line)
                if time is not None:
                    time_line = old_line

                maps = re.match(r'CheckOutputMaps=', old_line)
                if maps is not None:
                    maps_line = old_line

        with open('Alteck' + str(ev) + '.run', 'rt') as file:
            filedata = file.read()
            if work_line is not None:
                filedata = filedata.replace(work_line, str("Work Directory=" + new_dir + '\n'))

            if map_line is not None:
                filedata = filedata.replace(map_line, str("Map Directory=" + new_dir + map_dir_end + '\n'))

            if rain_line is not None:
                filedata = filedata.replace(rain_line, str("Rainfall Directory=" + new_dir + rain_dir_end + '\n'))

            if rain_f_line is not None:
                filedata = filedata.replace(rain_f_line, str("Rainfall file=Event_" + str(ev) + ".txt" + '\n'))

            if res_line is not None:
                filedata = filedata.replace(res_line, str("Result Directory=" + new_dir + res_dir_end + '\n'))

            if end_time_line is not None:
                filedata = filedata.replace(end_time_line, str("End time=" + str(end_time[events.index(ev)]) + '\n'))

            if time_line is not None:
                filedata = filedata.replace(time_line, str("Timestep=" + str(ts_new[events.index(ev)]) + '\n'))

            if maps_line is not None:
                filedata = filedata.replace(maps_line, str("CheckOutputMaps=" + map_zeroes + '\n'))

        with open('Alteck' + str(ev) + '.run', 'w') as nfile:
            nfile.write(filedata)


if __name__ == "__main__":
    update_runfile_directories(lisem_runs)

