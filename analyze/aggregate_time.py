import os

dataset=["email-Enron-full", "email-Eu-full", 
        "contact-high-school", "contact-primary-school",
        "NDC-classes-full", "NDC-substances-full",
        "tags-ask-ubuntu", "tags-math-sx", "threads-ask-ubuntu", "coauth-MAG-Geology-full", "coauth-MAG-History-full"]

alphaset1=["0.0000", "0.2500", "0.5000", "1.0000", "2.0000"]
alphaset=["0.00", "0.25", "0.50", "1.00", "2.00"]
betaset=["-1.00", "-0.50", "-0.25", "0.00", "0.25", "0.50", "1.00"]
ri = 5

def read_time(inputdirlist, conditionflag=False):
    samplingtime = 0
    propertytime = 0
    sizeavgtime = 0

    for inputdir in inputdirlist:
        with open(inputdir + "time.txt", "r") as f:
            for line in f.readlines():
                samplingtime += int(line.rstrip().split(" ")[0])
        if conditionflag:
            assert os.path.isfile(inputdir + "time_property.txt")
        if os.path.isfile(inputdir + "time_property.txt"):
            with open(inputdir + "time_property.txt", "r") as f:
                for line in f.readlines():
                    propertytime += int(line.rstrip().split(" ")[0])
        if conditionflag:
            assert os.path.isfile(inputdir + "time_size_avg.txt")
        if os.path.isfile(inputdir + "time_size_avg.txt"):
            with open(inputdir + "time_size_avg.txt", "r") as f:
                for line in f.readlines():
                    sizeavgtime += int(line.rstrip().split(" ")[0])

    return samplingtime, propertytime, sizeavgtime

for data in dataset:
    # baseline
    for algo in ["ns/add_global_deg_0.0000", "ns/add_global_deg_1.0000", "rw/rw_c_1.00_0.30", "ns/add_ordered", "es/add_global_deg_min_0.0000", "ff/ff_c_0.51_0.20", "tihs"]:
        if data == "coauth-MAG-Geology-full" and algo == "ff/ff_c_0.51_0.20":
            inputdir = "../results/{}/{}/{}/".format(algo, data, 6)
        else:
            inputdir = "../results/{}/{}/{}/".format(algo, data, ri)
        outputdir = "../results_time/{}/{}/".format(algo, data)
        if os.path.isdir(outputdir) is False:
            os.makedirs(outputdir)
        with open(inputdir + "time_sampling_portion.txt", "r") as f:
            lines = []
            for line in f.readlines():
                lines.append(line.rstrip().split(" ")[0])
            assert len(lines) == 10
        time1 = int(lines[8]) # 90%
        with open(outputdir + "samplingtime.txt", "w") as f:
            f.write(str(time1) + "\n")

    # midas
    inputdirlist = ["../results/es/add_global_deg_min_{}/{}/{}/".format(alpha, data, ri) for alpha in alphaset1]
    outputdir = "../results_time/midas/{}/".format(data)
    if os.path.isdir(outputdir) is False:
        os.makedirs(outputdir)
    time1, time2, time3 = read_time(inputdirlist, conditionflag=True)
    with open(outputdir + "samplingtime.txt", "w") as f:
        f.write(str(time1) + "\n")
    with open(outputdir + "propertytime.txt", "w") as f:
        f.write(str(time2 + time3) + "\n")
    with open(outputdir + "totaltime.txt", "w") as f:
        f.write(str(time1 + time2 + time3) + "\n")

    # midasB
    inputdirlist = []
    for alpha in alphaset:
        for beta in betaset:
            inputdirlist.append("../results/essz/add_global_deg_min_{}_{}/{}/{}/".format(alpha, beta, data, ri))
    outputdir = "../results_time/midasB/{}/".format(data)
    if os.path.isdir(outputdir) is False:
        os.makedirs(outputdir)
    time1, time2, time3 = read_time(inputdirlist, conditionflag=True)
    with open(outputdir + "samplingtime.txt", "w") as f:
        f.write(str(time1) + "\n")
    with open(outputdir + "propertytime.txt", "w") as f:
        f.write(str(time2 + time3) + "\n")
    with open(outputdir + "totaltime.txt", "w") as f:
        f.write(str(time1 + time2 + time3) + "\n")

