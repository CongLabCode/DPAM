import argparse
import docker
import os,sys

def check_image_exists(image_name):
    client = docker.from_env()
    try:
        client.images.get(image_name)
        return True
    except docker.errors.ImageNotFound:
        return False


def check_databases(databases_dir):
    path = os.getcwd()
    flag = 1
    if not os.path.exists(databases_dir):
        print (databases_dir, 'does not exist')
        flag = 0
    else:
        missing = []
        with open(f'{databases_dir}/all_files') as f:
            all_files = f.readlines()
        all_files = [i.strip() for i in all_files]
        for fn in all_files:
            if not os.path.exists(f'{databases_dir}/{fn}'):
                missing.append(fn)
        if missing:
            flag = 0
            with open('dpam_databases_missing_files','w') as f:
                f.write('\n'.join(missing)+'\n')
            print(f"Files missing for databases. Please check {path}/dpam_databases_missing_files for details")
        else:
            if os.path.exists('dpam_databases_missing_files'):
                os.system('rm dpam_databases_missing_files')
    return flag

def check_inputs(input_dir,dataset):
    flag = 1
    if not os.path.exists(input_dir):
        flag = 0
        print('Error!', input_dir, 'does not exist.')
    else:
        if os.path.exists(f'{input_dir}/{dataset}') and os.path.exists(f'{input_dir}/{dataset}_struc.list'):
            with open(f'{input_dir}/{dataset}_struc.list') as f:
                alist = f.readlines()
            alist = [i.strip() for i in alist]
            missing = []
            for name in alist:
                if not os.path.exists(f'{input_dir}/{dataset}/{name}.cif') and not os.path.exists(f'{input_dir}/{dataset}/{name}.pdb'):
                    missing.append([dataset, name, ':PDB/CIF missing'])
                if not os.path.exists(f'{input_dir}/{dataset}/{name}.json'):
                    missing.append([dataset, name, ':PAE json missing'])
            if missing:
                flag = 0
                with open(f'{input_dir}/dpam_{dataset}_inputs_missing_files','w') as f:
                    for i in missing:
                        f.write(' '.join(i)+'\n')
                print(f'Error! Please check {input_dir}/dpam_{dataset}_inputs_missing_files for details')
        if not os.path.exists(f'{input_dir}/{dataset}'):
            flag = 0
            print('Error!', dataset, 'containing PDB/CIF and PAE does not exist.')
        if not os.path.exists(f'{input_dir}/{dataset}_struc.list'):
            flag = 0
            print('Error!', f'{input_dir}/{dataset}_struc.list for targets does not exist.')
        return flag

         

def run_docker_container(image_name, databases_dir, input_dir, dataset, threads, log_file):
    client = docker.from_env()
    wdir = f'/home/'+input_dir.split('/')[-1]

    # Mount the directories to the container
    volumes = {
        databases_dir: {'bind': '/mnt/databases', 'mode': 'ro'},
        input_dir: {'bind': wdir, 'mode': 'rw'}
    }

    container = client.containers.run(image_name, detach=True, volumes=volumes, working_dir=wdir, command='tail -f /dev/null')

    # Example of running a script inside the container
    # Modify as needed for your specific script execution
    try:
        exec_log = container.exec_run(f"/bin/bash -c 'run_dpam.py {dataset} {threads}'", stdout=False, stderr=True)
        final_status = f'DPAM run for {dataset} under {input_dir} done\n'
    except:
        final_status = f'DPAM run for {dataset} under {input_dir} failed\n'

    with open(log_file, 'w') as file:
        file.write(exec_log.output.decode())
        file.write(final_status)
    # Stop the container after the script execution
    container.stop()

    # Optionally, remove the container if not needed anymore
    container.remove()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run a DPAM docker container.")
    parser.add_argument("--databases_dir", help="Path to the databases directory to mount (required)", required=True)
    parser.add_argument("--input_dir", help="Path to the input directory to mount (required)", required=True)
    parser.add_argument("--dataset", help="Name of dataset (required)", required=True)
    parser.add_argument("--image_name", help="Image name", default="dpam")
    parser.add_argument("--threads", type=int, default=os.cpu_count(), help="Number of threads. Default is to use all CPUs")
    parser.add_argument("--log_file", help="File to save the logs")

    args = parser.parse_args()

    image_flag = check_image_exists(args.image_name)
    if not image_flag:
        print(args.image_name, 'does not exist!')
        sys.exit(1)

    db_flag = check_databases(args.databases_dir)
    if db_flag == 0:
        print("Databases are not complete")
        sys.exit(1)

    input_flag = check_inputs(args.input_dir,args.dataset)
    if input_flag == 0:
        print('Error(s)! Inputs missing')
        sys.exit(1)

    if '/' != args.input_dir[0]:
        path = os.path.join(os.getcwd(), args.input_dir) 
        input_dir = os.path.abspath(path)
    else:
        input_dir = os.path.abspath(args.input_dir)

    if '/' != args.databases_dir[0]:
        path = os.path.join(os.getcwd(), args.databases_dir)
        databases_dir = os.path.abspath(path)
    else:
        databases_dir = os.path.abspath(args.databases_dir)
    
    if args.log_file is None: 
        log_file = input_dir + '/' + args.dataset + '_docker.log'
    else:
        log_file = args.log_file
    
    run_docker_container(args.image_name,databases_dir, input_dir, args.dataset, args.threads,log_file)
