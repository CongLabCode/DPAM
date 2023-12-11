#!/opt/conda/bin/python
import os, sys
import numpy as np
import tensorflow as tf
from tensorflow.python.client import device_lib

dataset = sys.argv[1]
local_device_protos = device_lib.list_local_devices()
gpus = [x.name for x in local_device_protos if x.device_type == 'GPU']
if not gpus:
    print("No GPUs found. Falling back to CPU.")
    config = tf.ConfigProto()
else:
    config = tf.ConfigProto(gpu_options = tf.GPUOptions(per_process_gpu_memory_fraction=0.9))
    config.gpu_options.allow_growth = True

fp = open('step16_' + dataset + '.list', 'r')
prots = []
for line in fp:
    words = line.split()
    prots.append(words[0])
fp.close()

all_cases = []
all_inputs = []
for prot in prots:
    if os.path.exists('step15/' + dataset + '/' + prot + '.data'):
        fp = open('step15/' + dataset + '/' + prot + '.data', 'r')
        for countl, line in enumerate(fp):
            if countl:
                words = line.split()
                all_cases.append([prot, words[0], words[1], words[2], words[3], words[17], words[18], words[19], words[20], words[21], words[22]])
                all_inputs.append([float(words[4]), float(words[5]), float(words[6]), float(words[7]), float(words[8]), float(words[9]), float(words[10]), float(words[11]), float(words[12]), float(words[13]), float(words[14]), float(words[15]), float(words[16])])
        fp.close()
total_case = len(all_cases)


def get_feed(batch_inputs):
    inputs = np.zeros((100, 13), dtype = np.float32)
    for i in range(100):
        for j, value in enumerate(batch_inputs[i]):
            inputs[i, j] = value
    feed_dict = {myinputs: inputs}
    return feed_dict

dense = tf.compat.v1.layers.dense
with tf.Graph().as_default():
    with tf.name_scope('input'):
        myinputs = tf.placeholder(dtype = tf.float32, shape = (100, 13))
    layers = [myinputs]
    layers.append(dense(layers[-1], 64, activation = tf.nn.relu))
    preds = dense(layers[-1], 2, activation = tf.nn.softmax)
    saver = tf.train.Saver()

    with tf.Session(config = config) as sess:
        saver.restore(sess, '/mnt/databases/domass_epo29')
        all_preds = []
        if total_case >= 100:
            batch_count = total_case // 100
            get_case = 0
            for i in range(batch_count):
                batch_inputs = all_inputs[i * 100 : i * 100 + 100]
                batch_preds = sess.run(preds, feed_dict = get_feed(batch_inputs))
                get_case += 100
                for j in range(100):
                    all_preds.append(batch_preds[j,:])
                if i % 1000 == 0:
                    print ('prediction for batch ' + str(i))

            remain_case = total_case - get_case
            add_case = 100 - remain_case
            batch_inputs = all_inputs[get_case:] + all_inputs[:add_case]
            batch_preds = sess.run(preds, feed_dict = get_feed(batch_inputs))
            for j in range(remain_case):
                all_preds.append(batch_preds[j,:])
        
        else:
            fold = 100 // total_case + 1
            pseudo_inputs = all_inputs * fold
            batch_inputs = pseudo_inputs[:100]
            batch_preds = sess.run(preds, feed_dict = get_feed(batch_inputs))
            for j in range(total_case):
                all_preds.append(batch_preds[j,:])

        prot2results = {}
        for prot in prots:
            prot2results[prot] = []
        for i in range(total_case):
            this_case = all_cases[i]
            this_input = all_inputs[i]
            this_pred = all_preds[i]
            prot = this_case[0]
            prot2results[prot].append([this_case[1], this_case[2], this_case[3], this_case[4], this_pred[1], this_input[3], this_input[4], this_input[5], this_input[6], this_input[7], this_input[8], this_input[9], this_input[10], this_input[11], this_input[12], this_case[5], this_case[6], this_case[7], this_case[8], this_case[9], this_case[10]])
        for prot in prots:
            if prot2results[prot]:
                rp = open('step16/' + dataset + '/' + prot + '.result', 'w')
                rp.write('Domain\tRange\tTgroup\tECOD_ref\tDPAM_prob\tHH_prob\tHH_cov\tHH_rank\tDALI_zscore\tDALI_qscore\tDALI_ztile\tDALI_qtile\tDALI_rank\tConsensus_diff\tConsensus_cov\tHH_hit\tDALI_hit\tDALI_rot1\tDALI_rot2\tDALI_rot3\tDALI_trans\n')
                for item in prot2results[prot]:
                    rp.write(f'{item[0]}\t{item[1]}\t{item[2]}\t{item[3]}\t{str(round(item[4], 4))}\t{str(item[5])}\t{str(item[6])}\t{str(item[7])}\t{str(item[8])}\t{str(item[9])}\t{str(item[10])}\t{str(item[11])}\t{str(item[12])}\t{str(item[13])}\t{str(item[14])}\t{item[15]}\t{item[16]}\t{item[17]}\t{item[18]}\t{item[19]}\t{item[20]}\n')
                rp.close()
            else:
                os.system('echo \'done\' > step16/' + dataset + '/' + prot + '.done')
