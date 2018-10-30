# -*- coding: utf-8 -*-
from model_g10_v2 import *
from run_manager import *
from SALib.sample import latin
from mlhs_v15_problem import *  # Defines the LHS sampling problem


def coupled_montecarlo_setup(smps=2, firstTimeStep=166):

    # Get and Save Latin Hypercube Matrix
    problem = get_problem()
    smps_matrix = latin.sample(problem, smps)
    # Normalized matrix (will be re-scaled by phd-model-process)
    saveLHSmatrix(smps_matrix)
    mc_upper = problem['upper']

    for sample_nr in range(1, smps+1):
        sample_vector = getInputVector(sample_nr-1, smps_matrix)
        run_coupled_sample(sample_nr, firstTimeStep,
                           sample_vector=sample_vector, mc_upper=mc_upper,
                           couple=True, montecarlo_couple=True)


def getInputVector(row, sample_matrix):
    """
    :param row: relevant sample row
    :param sample_matrix: numpy sample matrix
    :return: a numpy row with required input parameters
    """
    smp_vector = sample_matrix[row]
    return smp_vector


def run_coupled_sample(sample_nr, firstTimeStep,  # 166 -> 14/03/2016
                       sample_vector=None, mc_upper=None,
                       couple=False, montecarlo_couple=False):
    """
    Runs ONE montecarlo sample. Set-up is intended for BEACH+LISEM coupling

    :param sample_nr:
    :param firstTimeStep: This is only the very first simulation timestep
    :param sample_vector:
    :param mc_upper:
    :param couple:
    :param montecarlo_couple:
    :return:
    """
    origin = os.getcwd()
    print(origin)

    """
    Input parameters to start BEACH
    """
    params = ['z3_factor',
              'cZ0Z1',
              'cZ',
              'c_adr',
              'k_g',
              'gamma01', 'gammaZ',
              'f_transp',
              'f_evap',
              'f_oc', 'k_oc',
              'beta_runoff',
              'dt_50_aged',
              'dt_50_ab',
              'dt_50_ref',
              'epsilon_iso',
              'beta_moisture'
              ]
    runs = 0
    correct = False
    end_days_remain = deepcopy(lisem_runs)
    for period in range(len(lisem_runs)):
        os.chdir(origin)
        if runs == 0:
            runs += 1
            firstTimeStep = firstTimeStep

            if couple:
                lastTimeStep = lisem_runs[period]
                if montecarlo_couple:
                    param_values = sample_vector
                    upper = mc_upper
            else:
                lastTimeStep = 286
                param_values = np.loadtxt("best_vector.txt")
                upper = np.ones(len(param_values)).tolist()
        else:
            runs += 1
            firstTimeStep = lisem_runs[period - 1] - 1
            lastTimeStep = lisem_runs[period]

        print("Starting BEACH")
        # print("run = " + str(runs))
        print("Sample nr. = " + str(sample_nr))
        print("period = " + str(period))
        myAlteck16 = BeachModel("maps\\static\\clone_nom.map",
                                params, param_values, upper, period,
                                lisem_runs, end_days_remain, event_type,
                                map_list, hydro_row, ts_new,
                                staticDT50=False, correction=correct,
                                coupled=couple, montecarlo_couple=montecarlo_couple, sample_nr=sample_nr)


        dynamicModel = DynamicFramework(myAlteck16,
                                        firstTimestep=firstTimeStep,
                                        lastTimeStep=lastTimeStep  # -> 183
                                        )  # an instance of the Dynamic Framework
        dynamicModel.run()
        end_days_remain.pop(0)
        correct = True  # Will tell BEACH to import LISEM outputs instead of initial conditions on next run
        print("Finished")
        print(datetime.today().strftime('%Y-%m-%d %HH:%MM'))
        if not couple:
            break
        else:
            print("Preparing LISEM")
            # Modify run file (appropriate time-steps for recording maps)
            # total: 510 minutes; timestep: 1248 sec/timestep ,
            #   break-day: 208min x 60sec = 12480 sec;
            #       output interval = 12480/1248 = 10.

            os.chdir(os.getcwd() + r"\\LISEM")
            path_run = r"..\\Alteck" + str(lisem_runs[period]) + ".run"
            command_line = r"lisem -b -r " + str(path_run)
            args = shlex.split(command_line)

            print("Starting LISEM")
            code = subprocess.call(args)

            if code == -1073741819:
                # end_run = False
                print("Code was good: " + str(code))
                print("Storing LISEM outlet text files in dir: opL_out")

                src = origin + "\\res\\" + str(period) + "\\hydro.txt"
                dst = origin + "\\visuals_opL\\Ls" + str(sample_nr) + "\\" + str(period) + "\\hydro.txt"
                shutil.copyfile(src, dst)
                src = origin + "\\res\\" + str(period) + "\\total.txt"
                dst = origin + "\\visuals_opL\\Ls" + str(sample_nr) + "\\" + str(period) + "\\total.txt"
                shutil.copyfile(src, dst)

                if runs == len(lisem_runs):  # Do final run!
                    firstTimeStep = lisem_runs[period] - 1
                    # end_run = True
                    runs += 1
                    print("Starting last BEACH run")
                    print("run = " + str(runs))
                    print("period = " + str(period + 1))
                    os.chdir(origin)
                    myAlteck16 = BeachModel("maps\\static\\clone_nom.map",
                                            params, param_values, upper, period + 1,
                                            lisem_runs, end_days_remain, event_type,
                                            map_list, hydro_row, ts_new,
                                            staticDT50=False, correction=correct,
                                            coupled=couple, montecarlo_couple=montecarlo_couple, sample_nr=sample_nr)
                    dynamicModel = DynamicFramework(myAlteck16,
                                                    firstTimestep=firstTimeStep,
                                                    lastTimeStep=286  # -> 183
                                                    )  # an instance of the Dynamic Framework
                    dynamicModel.run()
                    # if len(end_days_remain) > 0:
                    # end_days_remain.pop(0)
                    break
                else:
                    continue
            else:
                print("Error code: " + str(code))
                print("Exiting coupled loop at smp nr. ", str(sample_nr))
                break

    print("All runs finished")
    print("Remaining days to couple: ", end_days_remain)
    print(datetime.today().strftime('%Y-%m-%d %HH:%MM'))


# Reduce "best" drainge coefficient
# correct_factor = .5 * .5 * .5
# new_values = []
# for i in range(len(names)):
#     if i != names.index("c_adr"):
#         new_values.append(best_values[i])
#     else:
#         new_values.append(best_values[i] * correct_factor)

# Inspect lisem inputs
#  aguila --scenarios='{0,1,2,3,4}' ksat1


if __name__ == "__main__":
    mc_samples = 2
    create_dirs(mc_samples)
    coupled_montecarlo_setup(smps=mc_samples, firstTimeStep=166)
