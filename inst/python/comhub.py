import warnings

warnings.filterwarnings('ignore')

"""
# =============================================================================
# ComHub
# Author: Julia Åkesson
#
# Copyright 2019 Julia Åkesson
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#  http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# =============================================================================
"""
import os
import numpy as np
import pandas as pd
import scipy.stats as sts
#import matplotlib.pyplot as plt
#import seaborn as sns

class comhub:
    """
    ComHub is a tool to make hub predictions. ComHub identifies hubs in GRNs by combining the results of a compendium of GRN predictions. ComHub selects an optimal threshold for the number of edges to include from the GRN predicitons. For each GRN prediction the outdegree of each regulator is calculated, before averaging the outdegree over all GRN predictions. The output is a list of regulators ranked on outdegree.

    Run:
    c = comhub(network_name, methods=['aracne', 'clr_R', 'pcc', 'elasticnet_bootstrap', 'tigress_R', 'genie3'], expression_data=None, transcription_factors=None, gold_standard=None)

    Run with MATLAB version of CLR and TIGRESS:
    c = comhub(network_name, methods=['aracne', 'clr', 'pcc', 'elasticnet_bootstrap', 'tigress', 'genie3'], expression_data=None, transcription_factors=None, gold_standard=None)
    """
    def __init__(self, network_name, methods=['aracne', 'clr_R', 'pcc', 'elasticnet_bootstrap', 'tigress_R', 'genie3'], expression_data=None, transcription_factors=None, gold_standard=None):

        self.methods = methods
        self.network_name = network_name

        if expression_data is None:
            self.expression_data = './data/'+network_name+'_expression_data.tsv'
        else:
            self.expression_data = expression_data

        if transcription_factors is None:
            self.transcription_factors = './data/'+network_name+'_transcription_factors.tsv'
        else:
            self.transcription_factors = transcription_factors

        if gold_standard is None and os.path.exists('./data/'+network_name+'_gold_standard.tsv'):
            self.gold_standard = './data/'+network_name+'_gold_standard.tsv'        
        elif gold_standard is not None:
            self.gold_standard = gold_standard        
           
        

    def run_methods(self, network_cutoff=100000, nstepsLARS=5, matlab=False, bootstrap=100, parallel=True):
        """
        Runs any combination of the network inference methods: aracne, clr_R, clr, pcc, elasticnet_bootstrap,
        tigress, tigress_R, and genie3. The methods should be specified when initiating comhub.

        run:
        c.run_methods(network_cutoff=100000)
        """
    
        snavel = self.pcc(network_cutoff=network_cutoff)
        
        return(snavel)

    def aracne(self, aracne_mat, network_cutoff=100000):
        """
        ARACNE “algorithm for the reconstruction of accurate cellular networks”
        1) Computes a mutual information score for each regulator-target interaction
        2) Uses data processing inequality to remove indirect interactions.

        ARACNE is implemented using the r-package minet.
        Additional dependencies:
        R, minet, rpy2

        Run:
        net = c.aracne(network_cutoff=100000)

        Reference:
        Margolin,A.A. et al. (2006) ARACNE: An algorithm for the reconstruction of gene regulatory networks in a mammalian cellular context. BMC Bioinformatics, 7(Suppl1), 1–15.
        """
        aracne_mat = pd.DataFrame(aracne_mat)
        tfs = pd.DataFrame(self.transcription_factors).astype(str)
        tfs = tfs.astype(str)
        tfs = tfs[tfs.isin(aracne_mat.columns)].dropna()
      

        aracne_mat = aracne_mat.loc[list(tfs.iloc[:, 0])]
        aracne_mat = aracne_mat.dropna()
        
        net = aracne_mat.reset_index().melt(id_vars=['index'])
        net.columns = ['TF', 'target', 'confidence']
        net_sort = net.sort_values(by='confidence', ascending=False).iloc[:network_cutoff, :]
        net_sort['confidence'] = pd.to_numeric(net_sort['confidence'])
        net_sort = net_sort[net_sort.confidence != 0]
        print('saving network')
        net_sort.to_csv('./networks/'+self.network_name+'/aracne_network.tsv', sep='\t', index=False, header=None)
        print('Done')
        return tfs
        #return net_sort

    def clr_R(self, network_cutoff=100000):
        """
        CLR “context likelihood of relatedness”
        1) Computes a mutual information score for each regulator-target interaction
        2) filters interactions not significantly above the “background” distribution of MI scores.

        CLR is implemented using the r-package minet.
        Additional dependencies:
        R, minet, rpy2

        Run:
        net = c.clr_R(network_cutoff=100000)

        Reference:
        Faith,J.J. et al. (2007) Large-scale mapping and validation of Escherichia coli transcriptional regulation from a compendium of expression profiles. PLoS Biol., 5, 0054–0066.
        """
        import rpy2.robjects as robjects
        from rpy2.robjects.packages import importr
        minet = importr('minet')
        
        ##run aracne in r
        readtable = robjects.r['read.table']
        exp = readtable(self.expression_data, header=True, sep='\t')
        print('building mutual information matrix')
        buildmim = robjects.r['build.mim']
        mim = buildmim(exp, estimator="spearman")
        print('running clr')
        clr = robjects.r['clr']
        clr_matrix = clr(mim)
        print('saving clr matrix')
        writetable = robjects.r['write.table']
        save_file = './networks/'+self.network_name+'/clr_matrix.csv'
        writetable(clr_matrix, file=save_file)
        
        ##Create network
        print('creating network')
        clr_matrix = pd.read_csv(save_file, sep=' ')
        exp = pd.read_csv(self.expression_data, sep='\t')
        clr_matrix.columns = exp.columns
        clr_matrix.index = exp.columns
        
        tfs = pd.read_csv(self.transcription_factors, header=None).astype(str)
        tfs = tfs[tfs.isin(clr_matrix.columns)].dropna()
        
        clr_matrix = clr_matrix.loc[list(tfs.iloc[:, 0])]
        clr_matrix = clr_matrix.dropna()
        
        net = clr_matrix.reset_index().melt(id_vars=['index'])
        net.columns = ['TF', 'target', 'confidence']
        net_sort = net.sort_values(by='confidence', ascending=False).iloc[:network_cutoff, :]
        net_sort['confidence'] = pd.to_numeric(net_sort['confidence'])
        net_sort = net_sort[net_sort.confidence != 0]
        print('saving network')
        net_sort.to_csv('./networks/'+self.network_name+'/clr_R_network.tsv', sep='\t', index=False, header=None)
        print('Done')
        return net_sort

    def clr(self, network_cutoff=100000, matlab=False):
        """
        CLR “context likelihood of relatedness”
        1) Computes a mutual information score for each regulator-target interaction
        2) filters interactions not significantly above the “background” distribution of MI scores.

        CLR calls MATLAB functions using either octave or MATLAB. MATLAB is recommended for faster performance.
        Additional dependencies:
        scikit-learn, Octave: octave, oct2py, MATLAB: matlab, matlab.engine

        Run:
        net = c.clr(network_cutoff=100000, matlab=False)

        Reference:
        Faith,J.J. et al. (2007) Large-scale mapping and validation of Escherichia coli transcriptional regulation from a compendium of expression profiles. PLoS Biol., 5, 0054–0066.
        """
        from sklearn.metrics import mutual_info_score

        def calc_mi(x, y, bins):
            c_xy = np.histogram2d(x, y, bins)[0]
            mi = mutual_info_score(None, None, contingency=c_xy)
            return mi

        #print('Reading files with expression data and transcription factors')
        exp = self.expression_data
        tf = self.transcription_factors

        #tf=pd.DataFrame(exp.columns)
        number_of_genes = len(exp.columns)
        exp_array = np.array(exp)

        print('calcualting mutual information')
        mi_full = np.zeros((number_of_genes, number_of_genes))
        for i in range(number_of_genes):
            for j in range(number_of_genes):
                mi_full[i, j] = calc_mi(x=exp_array[:, i], y=exp_array[:, j], bins=10)

        mi_matrix = mi_full
        genes = np.array(exp.columns)
        pd.DataFrame(mi_matrix).to_csv('./networks/'+self.network_name+'/clr_matrix.csv')

        ##Matlab part#
        if matlab:
            import matlab.engine
            print('running clr in matlab')
            eng = matlab.engine.start_matlab()
            eng.addpath('./bin/')
            m = matlab.double(mi_matrix.tolist())
            clr_matrix = eng.clr_octave(m, 'normal')
            eng.quit()
        else:
            from oct2py import octave
            print('running clr in octave')
            octave.addpath('/home/dirk/commifier/bin/')
            clr_matrix = octave.clr_octave(mi_matrix, 'normal')

        ## postive values should mean that MI is significantly above background
        ## and negative values should mean that values are significantly below background.

        clr = pd.DataFrame(clr_matrix, columns=genes, index=genes)

        tfs = np.array(tf).flatten()
        tfs = [str(x) for x in tfs]

        clr_dir = clr.loc[tfs, :] #Only keep rows that are TFs to get directed network, see above.

        net = clr_dir.reset_index().melt(id_vars=['index'])
        net.columns = ['TF', 'target', 'confidence']
        net.confidence = [float(x) for x in list(net.confidence)]
        net = net[net.TF != net.target]
        net_sort = net.sort_values(by='confidence', ascending=False).iloc[:network_cutoff, :]

        print('saving network')
        net_sort.to_csv('./networks/'+self.network_name+'/clr_network.tsv', index=False, header=None, sep='\t')
        print('Done')
        return net_sort

    def pcc(self, network_cutoff=100000):
        """
        Absolute value of the Pearson correlation coefficient (PCC)
        Regulator-target interactions are ranked based on the absolute value of the PCC.

        Run:
        net = c.pcc(network_cutoff=100000)
        """
        exp = self.expression_data.T
        tfs = self.transcription_factors
        tfs = np.array(tfs).flatten().astype(str)
        tfs = [tf for tf in tfs if tf in exp.columns]
        tfs_index = np.where(np.isin(exp.columns,tfs) == True)[0]
        tfs_new_order = exp.columns[tfs_index]
        
        print('calculating pcc')
        corrmat = np.corrcoef(np.array(exp.T))
        corrmat = np.abs(corrmat[tfs_index,:])
        corrmat = pd.DataFrame(corrmat, index=tfs_new_order, columns=exp.columns)
        
        print('making network')
        net = corrmat.reset_index().melt(id_vars=['index'])
        net.columns = ['TF', 'target', 'confidence']
        net = net[net.TF != net.target]
        net_sort = net.sort_values(by='confidence', ascending=False).iloc[:network_cutoff, :]
        
        return net_sort

    def elasticnet_bootstrap(self, bootstrap=100, parallel=True, network_cutoff=100000,name=''):
        """
        Bootstrap Elastic Net

        Additional dependencies:
        scikit-learn, joblib (parallel)

        Run:
        net = c.elasticnet_bootstrap(bootstrap=100, parallel=True, network_cutoff=100000)

        Reference:
        Zou,H. and Hastie,T. (2005) Regularization and variable selection via the elastic-net. J. R. Stat. Soc., 67, 301–320.
        """
        import sklearn.linear_model as lm
        from sklearn.utils import resample
        if parallel:
            from joblib import Parallel, delayed

        def run_en(exp, tfs_index):
            targets = np.transpose(resample(np.array(exp.T)))
            tfs = np.transpose(targets[tfs_index])
            model = lm.ElasticNetCV(n_jobs=-1, cv=3)
            coef_mat = []
            for i, j in enumerate(targets):
                #print(str(i) + ' of ' + str(len(targets.T)))
                model.fit(tfs, targets[i])
                coef_mat.append(model.coef_)
            coef_mat = (np.array(coef_mat) > 0)*1
            return coef_mat

        #read data
        print('reading data')
        exp = self.expression_data
        tfs = self.transcription_factors
        print(exp.shape)
        tfs = np.array(tfs).flatten().astype(str)
        tfs = [tf for tf in tfs if tf in exp.index]
        print(tfs)
        tfs_index = np.where(np.in1d(np.array(exp.index), np.array(tfs)))[0]
        tfs_new_order = exp.index[tfs_index]

        #ElasticNet with bootstrap
        print('running Elastic Net')

        if parallel:
            result = Parallel(n_jobs=-1)(delayed(run_en)(exp, tfs_index) for i in range(bootstrap))
            result_sum = np.sum(result, axis=0)
        else:
            result_sum = np.zeros((len(exp), len(tfs_index)))
            for i in range(bootstrap):
                result = run_en(exp, tfs_index)
                result_sum = result_sum + result
       
        coef_mat = pd.DataFrame(result_sum/bootstrap, columns=tfs_new_order, index=exp.index).T
       
        print('Making network')
        net = coef_mat.reset_index().melt(id_vars=['index'])
        net.columns = ['TF', 'target', 'confidence']
        net = net[net.TF != net.target]
        net_sort = net.sort_values(by='confidence', ascending=False).iloc[:network_cutoff, :]
        net_sort = net_sort[net_sort.confidence != 0]

        #net_sort.to_csv('./networks/'+self.network_name+'/elasticnet_bootstrap_network.tsv', index=False, header=None, sep='\t')
        print('Done')
        return net_sort

    def tigress_R(self, network_cutoff=100000, nstepsLARS=5):
        """
        TIGRESS “trustful inference of gene regulation with stability selection”
        1) Least angle regression (LARS)
        2) stability selection

        TIGRESS is implemented using the r-package tigress (Downloaded from https://github.com/jpvert/tigress).
        Additional dependencies:
        R, tigress, rpy2

        Run:
        net = tigress_R(network_cutoff=100000, nstepsLARS=5)

        Reference:
        Haury, A.C. et al. (2012) TIGRESS: Trustful Inference of Gene REgulation using Stability Selection. BMC Syst. Biol., 6, 1–17.
        """
        
        import rpy2
        from rpy2.robjects import r, pandas2ri
        import rpy2.robjects as robjects
        from rpy2.robjects.packages import importr
        tigress = importr('tigress')
        
        print('running tigress')
        readtable = robjects.r['read.table']
        exp = readtable(self.expression_data, header=True, sep='\t', check_names=False)
        exp_targets = pd.read_csv(self.expression_data, sep='\t')
        tfs = pd.read_csv(self.transcription_factors, sep='\t', header=None)
        tfs = np.array(tfs).flatten().astype(str)
        tfs = [tf for tf in tfs if tf in exp_targets.columns]
        tfs = rpy2.robjects.vectors.StrVector(tfs)
        targets = rpy2.robjects.vectors.StrVector(exp_targets.columns)

        tigress = robjects.r['tigress']
        edgepred = tigress(exp, tflist=tfs, targetlist=targets, nstepsLARS=nstepsLARS)
        edgepred = edgepred[nstepsLARS-1]
        writetable = robjects.r['write.table']
        print('saving tigress matrix')
        save_file = './networks/'+self.network_name+'/tigress_matrix.csv'
        writetable(edgepred, file=save_file)

        ##Create network
        print('creating network')
        tigress_mat = pd.read_csv(save_file, sep=' ')
        #tigress_mat.columns = tigress_mat.columns.str.replace('X','')
        
        net = tigress_mat.reset_index().melt(id_vars=['index'])
        net.columns = ['TF', 'target', 'confidence']
        net_sort = net.sort_values(by='confidence', ascending=False).iloc[:network_cutoff, :]
        net_sort['confidence'] = pd.to_numeric(net_sort['confidence'])
        net_sort = net_sort[net_sort.confidence != 0]
        print('saving network')
        net_sort.to_csv('./networks/'+self.network_name+'/tigress_R_network.tsv', sep='\t', index=False, header=None)
        print('Done')
        return net_sort

    def tigress(self, network_cutoff=100000, matlab=False):
        """
        TIGRESS “trustful inference of gene regulation with stability selection”
        1) Least angle regression (LARS)
        2) stability selection

        TIGRESS calls MATLAB functions using either Octave or MATLAB. MATLAB is recommended for faster performance.
        Additional dependencies:
        Octave: octave, oct2py, MATLAB: matlab, matlab.engine

        Run:
        net = tigress(network_cutoff=100000, matlab=False)

        Reference:
        Haury, A.C. et al. (2012) TIGRESS: Trustful Inference of Gene REgulation using Stability Selection. BMC Syst. Biol., 6, 1–17.
        """
        if matlab:
            import matlab.engine
            print('running tigress in matlab')
            path = os.getcwd()
            eng = matlab.engine.start_matlab()
            eng.addpath('./bin/tigress/')
            net = eng.tigress_matlab(path, self.network_name, self.expression_data, self.transcription_factors, str(network_cutoff))
            eng.quit()
            print('Done')
        else:
            from oct2py import octave
            print('running tigress in octave')
            path = os.getcwd()
            octave.addpath('./bin/tigress/')
            net = octave.tigress_octave(path, self.network_name, self.expression_data, self.transcription_factors, str(network_cutoff))
            print('Done')
        return net

    def genie3(self, network_cutoff=100000):
        """
        GENIE3 “gene network inference with ensemble of trees”
        Decomposes the network inference into different feature selection problems and applies tree-based ensemble methods on each sub-problem.

        Run:
        net = c.genie3(network_cutoff=100000)

        Reference:
        Huynh-Thu,V.A. et al. (2010) Inferring regulatory networks from expression data using tree-based methods. PLoS One, 5, 1–10.
        """
        from bin.GENIE3 import loadtxt
        from bin.GENIE3 import GENIE3
        from bin.GENIE3 import get_link_list

        print('reading data')
        data = loadtxt(self.expression_data, skiprows=1)
        f = open(self.expression_data)
        genes = f.readline()
        f.close()
        genes = genes.rstrip('\n').split('\t')
        tf = pd.read_csv(self.transcription_factors, header=None)
        tf = [str(x) for x in list(np.array(tf).flatten())]

        print('running GENIE3')
        vim = GENIE3(data, gene_names=genes, regulators=tf)
        print('saving network')
        if not network_cutoff:
            net = get_link_list(vim, gene_names=genes, regulators=tf, maxcount='all',
                                file_name='./networks/'+self.network_name+'/genie3_network.tsv')
        else:
            net = get_link_list(vim, gene_names=genes, regulators=tf, maxcount=network_cutoff,
                                file_name='./networks/'+self.network_name+'/genie3_network.tsv')
        print('Done')
        return net

    def get_tf_outdegree(self, network_files, names, edge_cutoff=100000):
        """
        Calculates the outdegree of regulators in a set of GRN predicitons for a specified edge threshold.
        Reads networks named './networks/{network_name}/{method}_network.tsv', if not a list of network files is specified.
        Run:
        tf = get_tf_outdegree(edge_cutoff=100000, network_files=None)
        """
        edge_cutoff = int(edge_cutoff)
        if not network_files:
            network_files = ['./networks/' + self.network_name + '/' + method + '_network.tsv' for method in self.methods]
        tf_outdegree_all = pd.DataFrame()
        for result, method in zip(network_files, names):
            method_network = result
            if len(method_network) > edge_cutoff:
                method_network = method_network.iloc[:edge_cutoff, :]
            tf_outdegree = pd.DataFrame(method_network.groupby(method_network.iloc[:, 0].tolist()).size().sort_values(ascending=False),
                                        columns=[method])
            tf_outdegree_all = tf_outdegree_all.join(tf_outdegree, how='outer')
        tf_outdegree_all = tf_outdegree_all.fillna(0).astype(int)
        return tf_outdegree_all

    def pairwise_correlation(self, network_files, names, edge_range=[500, 1000, 2000, 3000, 4000, 5000, 7000, 10000, 15000, 20000, 50000, 80000, 100000], plot=True, fig_name=''):
        """
        Identifies an optimal edge threshold by assessing the pairwise correlation
        among GRN predictions for a range of edge thresholds.
        Outputs a figure named "comhub/results/{network_name}/pairwise_correlation.png".
        Run:
        edge_cutoff = pairwise_correlation(self, edge_range=[500, 1000, 2000, 3000, 4000, 5000, 7000, 10000, 15000, 20000, 50000, 80000, 100000], plot=True, fig_name='')
        """
        pair_pcc = pd.DataFrame()
        for nb_edges in edge_range:
            tf_outdegree = self.get_tf_outdegree(network_files = network_files, names = names, edge_cutoff=nb_edges)
            #Pairwise correlation.
            pcc_all = pd.DataFrame(columns=[str(nb_edges)])
            pcc_pval_all = pd.DataFrame(columns=[str(nb_edges)])
            for i in range(len(tf_outdegree.columns)):
                for j in range(len(tf_outdegree.columns)-i-1):
                    pcc = sts.pearsonr(tf_outdegree.iloc[:, i], tf_outdegree.iloc[:, i+j+1])
                    pcc_all = pcc_all.append(pd.DataFrame(pcc[0],
                                                          columns=[tf_outdegree.columns[i]+'_'+tf_outdegree.columns[i+j+1]],
                                                          index=[str(nb_edges)]).T)
                    pcc_pval_all = pcc_pval_all.append(pd.DataFrame(pcc[1],
                                                                    columns=[tf_outdegree.columns[i]+'_'+tf_outdegree.columns[i+j+1]],
                                                                    index=[str(nb_edges)]).T)
            pair_pcc = pair_pcc.join(pcc_all, how='outer')

        edge_cutoff = pair_pcc.mean().idxmax()

        #if plot:
            # pair_pcc.columns = [int(x) for x in pair_pcc.columns]
            # pair_pcc = pair_pcc.melt().dropna()
            # pair_pcc.columns = ['interactions', 'pcc']
            # pair_pcc.pcc = pair_pcc.pcc.astype(float)
            # 
            # color_pallet = "bright"
            # plt.style.use('seaborn-ticks')
            # sns.set_color_codes(color_pallet)
            # fig, ax = plt.subplots(ncols=1, nrows=1, facecolor='w', edgecolor='k', figsize=(15, 15), dpi=75)
            # plt.rc('font', size=36)
            # plt.rc('ytick', labelsize=32)
            # plt.rc('xtick', labelsize=32)
            # 
            # sns.lineplot(x='interactions', y='pcc', data=pair_pcc, ci=95, marker='D', markersize=7)
            # ax.axvline(x=int(edge_cutoff), c='k', linestyle='--', linewidth=2)
            # ax.set_xscale('log')
            # sns.despine()
            # fig.savefig('./results/'+self.network_name+'/pairwise_correlation'+fig_name+'.png', dpi=200, bbox_inches="tight")

        return int(edge_cutoff)

    def community(self, tf_outdegree, save_csv=True, output_name=''):
        """
        The outdegree of each regulator is averaged over the method predictions.
        Run:
        community = community(tf_outdegree, save_csv=True, output_name='')
        """
        community = np.mean(tf_outdegree, axis=1).sort_values(ascending=False)
        if save_csv:
            community.to_csv('./results/'+self.network_name+'/community'+output_name+'.tsv', sep='\t', header=False)
        return community

    def hubs(self, community, percentage_interactions=0.1, save_csv=True, output_name=''):
        """
        Identifies hubs among regulators. Top-ranked regulators standing
        for a certain percentage of the interactions in the network are identified as hubs.
        Run:
        hubs = hubs(community, percentage_interactions=0.1, save_csv=True, output_name='')
        """
        hubs = community.cumsum()[community.cumsum() < community.sum()*percentage_interactions]
        #hubs = community.index[:int(len(community)*percentage)]
        if save_csv:
            hubs.to_csv('./results/'+self.network_name+'/hubs'+output_name+'.tsv', sep='\t', header=False)
        return hubs

    def community_performance(self, community):
        """
        Evaluates the peformance of the community if a gold standard is available.
        The performance is evaluated with the Pearson correlation coefficient.
        Run:
        pcc, pval = community_performance(community)
        """
        gs = pd.read_csv('./data/'+self.network_name+'_gold_standard.tsv', header=None, sep='\t')
        gs_outdegree = pd.Series(gs.groupby(gs.loc[:, 0].tolist()).size().sort_values(ascending=False))
        community = community.loc[gs_outdegree.index]
        pcc, pval = sts.pearsonr(gs_outdegree, community)
        return pcc, pval

    def method_performance(self, edge_cutoff, plot=True):
        """
        Evaluates the performance of ComHub and each of the method predicitons.
        The performance is evaluated with the Pearson correlation coefficient.
        Outputs a figure named "comhub/results/{network_name}/method_performance_{edge_cutoff}.png"
        Run:
        performance = method_performance(edge_cutoff, plot=True)
        """
        def pearson_corr(x, Ymat):
            pcc = []
            pval = []
            for i in Ymat:
                pcof = sts.pearsonr(x, Ymat.loc[:, i])
                pcc.append(pcof[0])
                pval.append(pcof[1])
            return pcc, pval

        gs = pd.read_csv('./data/'+self.network_name+'_gold_standard.tsv', header=None, sep='\t')
        gs_outdegree = pd.Series(gs.groupby(gs.loc[:, 0].tolist()).size().sort_values(ascending=False))
        tf_outdegree = self.get_tf_outdegree(edge_cutoff=edge_cutoff)
        tf_outdegree = tf_outdegree.reindex(gs_outdegree.index)
        if tf_outdegree.isnull().sum().sum() > 0:
            tf_outdegree = tf_outdegree.dropna()
            print('Warning: all gold standard TFs are not present among the possible TFs')
        methods_pcc, methods_pval = pearson_corr(x=gs_outdegree, Ymat=tf_outdegree)

        community = self.community(tf_outdegree=tf_outdegree)
        community_pcc, community_pval = self.community_performance(community)
        pcc = methods_pcc + [community_pcc]
        pval = methods_pval + [community_pval]
        performance = pd.DataFrame([pcc, pval], columns=self.methods+['ComHub'], index=['pcc', 'pval']).T

        #if plot:
            # color_pallet = "bright"
            # plt.style.use('seaborn-ticks')
            # sns.set_color_codes(color_pallet)
            # 
            # fig, ax = plt.subplots(facecolor='w', edgecolor='k', figsize=(10, 8), dpi=75)
            # plt.rc('ytick', labelsize=32)
            # plt.rc('xtick', labelsize=32)
            # 
            # ax.bar(np.arange(0, len(pcc)), height=np.array(pcc), width=0.9)
            # ax.set_ylabel('PCC', fontsize=32)
            # 
            # plt.setp(ax, xticks=np.arange(len(self.methods+['ComHub'])), xticklabels=self.methods+['ComHub'])
            # sns.despine(offset=0, trim=False, bottom=True)
            # ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=True)
            # plt.xticks(rotation=90)
            # fig.savefig('./results/'+self.network_name+'/method_performance_'+str(edge_cutoff)+'.png', dpi=200, bbox_inches="tight")

        return performance

    def method_performance_edge_range(self, edge_range = [500, 1000, 2000, 3000, 4000, 5000, 7000, 10000, 15000, 20000, 50000, 80000, 100000], plot=True):
        def pearson_corr(x, Ymat):
            pcc = []
            pval = []
            for i in Ymat:
                pcof = sts.pearsonr(x, Ymat.loc[x.index, i].fillna(0))
                pcc.append(pcof[0])
                pval.append(pcof[1])
            return pcc, pval
        
        gs = pd.read_csv('./data/'+self.network_name+'_gold_standard.tsv', header=None, sep='\t')
        gs_outdegree = pd.Series(gs.groupby(gs.loc[:, 0].tolist()).size().sort_values(ascending=False))
                
        gs_pcc = pd.DataFrame()
        for nb_edges in edge_range:
            tf_outdegree = self.get_tf_outdegree(edge_cutoff=nb_edges)
            pcc = pd.DataFrame(pearson_corr(x=gs_outdegree, Ymat=tf_outdegree)[0])
            pcc.columns = [str(nb_edges)]
            gs_pcc = gs_pcc.join(pcc, how='outer')
        
        #if plot:
            # gs_pcc.columns = [int(x) for x in gs_pcc.columns]
            # gs_pcc = gs_pcc.melt()
            # gs_pcc.columns = ['Interactions', 'PCC']
            # 
            # color_pallet = "bright"
            # plt.style.use('seaborn-ticks')
            # sns.set_color_codes(color_pallet)
            # fig, ax = plt.subplots(ncols=1, nrows=1, facecolor='w', edgecolor='k', figsize=(15,15), dpi=75)
            # plt.rc('font', size=36)
            # plt.rc('ytick', labelsize=32)
            # plt.rc('xtick', labelsize=32)
            # 
            # sns.lineplot(x='Interactions', y='PCC', data=gs_pcc, ci=95, marker='D', markersize=7)
            # ax.set_xscale('log')
            # sns.despine()
            # fig.savefig('./results/'+self.network_name+'/goldstandard_correlation.png', dpi=200, bbox_inches="tight")
        return
