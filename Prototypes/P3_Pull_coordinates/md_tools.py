def mdp_create(file_name, new_parameters = None):      # for the .mdp-file
    if new_parameters is None:
        new_parameters = {}
    keys, parameters = read_mdp("base.mdp")
    orig_no = len(keys)
    for key in new_parameters.keys():
        if key not in parameters:
            keys.append(key)
        parameters[key] = new_parameters[key]
    file = open(file_name,"w")
    for par_count in range(len(keys)):
        if par_count == orig_no:
            file.write("\nCustom_parameters\n")
        key = keys[par_count]
        value = parameters[key]
        file.write("{:25} = {}\n".format(key,value))


def read_mdp(file_name):
    file =  open(file_name)
    keys=[]
    parameters = {}
    for row in file.read().split("\n"):
        if row == "":
            continue
        row = row.split("=")
        key = row[0].strip()
        value = row[1].strip()
        keys.append(key)
        parameters[key]=value
    file.close()
    return keys, parameters

if __name__ == "__main__":
    mdp_create("simple.mdp")
