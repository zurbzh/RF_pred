import sklearn
print sklearn.__version__  #it would be great if you had version 0.18

import numpy as np
import matplotlib.pyplot as plt
import csv

from numpy import loadtxt
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score

#the normal random shuffle split
def TT_split(features, labels, prop=0.9):
    nr_samples = features.shape[0]
    nr_train = int(nr_samples*prop)
    indx = range(nr_samples)
    np.random.shuffle(indx)
    features_temp = features[indx,:]
    labels_temp = labels[indx]

    train_X = features_temp[:nr_train,:]
    train_Y = labels_temp[:nr_train]
    test_X = features_temp[nr_train:,:]
    test_Y = labels_temp[nr_train:]
    return train_X,train_Y,test_X,test_Y

#split by giving an array of experiment numbers and the nr to be left to the test set
def project_split(features, labels, projects, experiment_nr, seqID):
    projects =np.array(projects)
    #print "test set size", np.shape(np.where(projects == experiment_nr)[0])
    assert(len(np.where(projects != experiment_nr)[0])>0)
    
    train_X = features[np.where(projects != experiment_nr)[0]]
    train_Y = labels[np.where(projects != experiment_nr)[0]]
    test_X = features[np.where(projects == experiment_nr)[0]]
    test_Y = labels[np.where(projects == experiment_nr)[0]]
    test_seqID = seqID[np.where(projects == experiment_nr)[0]]
    print np.shape(test_seqID)
    #print "PROJECT SPLIT: total,train, test sizes", np.shape(labels),np.shape(train_Y),np.shape(test_Y)
    return train_X,train_Y,test_X,test_Y,test_seqID


def column(matrix, i):
    return [row[i] for row in matrix]

#helper to get the project ID from the sequenceID (which is a string)
def seqID_to_projects(seqID):
    seq_projects = np.array([x.split("rr", 1) for x in seqID])
    projects = np.array(column(seq_projects,1))
    
    #print np.unique(projects), np.shape(np.unique(projects))
    exp_names = np.unique(projects)
    print "EXPERIMENT NAMES: "
    pr_nr = np.zeros(np.shape(projects))
    for i,name in enumerate(exp_names):
        print i, name, np.shape(np.where(projects==name))
        pr_nr[np.where(projects==name)[0]] = i
    print "\n\n\n"
    return pr_nr,exp_names
    
    
f=["data/DNA_1/concat_al2"]

#FOR AREA UNDER CURVE
roc_probas=[]
roc_labels=[]

#do you want to get precision-recal infor for each experiment?
print_out = True

for i,folder in enumerate(f):
    print "---------#################  LEAVE ONE OUT  ###################--------\n",folder
    if "concat_al2" in folder:
        seqID1 = loadtxt(folder[:-10]+"forward_al2/RCSU_id.txt",dtype=str)
        features1 = loadtxt(folder[:-10]+"forward_al2/matrix.txt")
        labels1 = loadtxt(folder[:-10]+"forward_al2/values.txt",dtype=int)
        
        features2 = loadtxt(folder[:-10]+"reverse_al2/matrix.txt")
        labels2 = loadtxt(folder[:-10]+"reverse_al2/values.txt",dtype=int)
        seqID2 = loadtxt(folder[:-10]+"reverse_al2/RCSU_id.txt",dtype=str)

        features = np.concatenate([features1,features2])
        labels = np.concatenate([labels1,labels2])
        seqID = np.concatenate([seqID1,seqID2])
 
        projects,exp_names=seqID_to_projects(seqID)
        
    else:
        print " NOT SUPPORTED!! We do LOO only on concat_al2 right now"
        assert(False)
        
    print "overall data shape: ",features.shape, labels.shape, "\n overall counts of classes ", np.bincount(labels)
    
    for LO in range(20):
        [train_X,train_Y,test_X,test_Y,test_seqID] = project_split(features, labels, projects,LO,seqID)

        cc = np.bincount(test_Y)
        if len(cc)==1: #if there is only one label to be counted
            print "\n\n LOO "+str(LO)+":  there are NO VIRUSES in the test set! Skipping this experiment."
            continue # end this iteration in the for loop

        #rf = RandomForestClassifier(class_weight="balanced")
        rf = RandomForestClassifier(n_estimators=100)
        #rf = RandomForestClassifier()
        rf.fit(train_X, train_Y)
        
        
        pred = rf.predict(test_X)
        prob = rf.predict_proba(test_X)
        
        xx = np.column_stack((test_seqID,prob, pred))
        
        with open("output.csv", "wb") as f:
    		writer = csv.writer(f)
    		writer.writerows(xx)
        
        
        roc_probas.append(prob)
        roc_labels.append(test_Y)     
        
        if print_out:
            print "\n\n LOO "+str(LO)+" train and test shapes after split ", train_X.shape, test_X.shape
            print "counts in test set", np.bincount(test_Y)
            report = classification_report(test_Y, pred, target_names=["not virus", "virus","not_classified"])
            print "LOO "+str(LO)+" report at 0.5: \n", report