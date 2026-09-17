### Figure 4d–f GSEApy
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gseapy as gp
import scanpy as sc
gp.__version__



labels=["DarkVsLight_Melanocyte"]
labels1=["Light"]
labels2=["Dark"]

for j in [0]:
    label=labels[j]
    label
    label1=labels1[j]
    label2=labels2[j]
    label
    label1
    label2

    # https://gseapy.readthedocs.io/en/master/singlecell_example.html
    adata = sc.read_h5ad(label + ".h5ad") # data from SeuratData::ifnb
    adata.obs.head()
    adata.X
    adata.layers['counts'] = adata.X # Save raw counts
    # preprocessing
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers['lognorm'] = adata.X
    adata.obs.groupby('orig.ident')['group'].value_counts()

    adata.obs['stim'] = pd.Categorical(adata.obs['group'], categories=[label2,label1], ordered=True)

    adata.obs.groupby('orig.ident')['stim'].value_counts()
    indices = adata.obs.sort_values(['stim']).index
    adata = adata[indices,:]

    # subset data
    bdata=adata.copy()

    # GSEA
    names = gp.get_library_name()
    print(names)
    #'MSigDB_Computational', 'MSigDB_Hallmark_2020'

    gene_sets=['GO_Biological_Process_2023','MSigDB_Hallmark_2020']
    setslabel=["GOBP","MSigHallmark"]

    for jj in [0,1]:
        set=gene_sets[jj]
        setlabel=setslabel[jj]
        set
        setlabel

        newpath = setlabel+'/'+label
        if not os.path.exists(newpath):
            os.makedirs(newpath)

        newpath = 'figures'
        if not os.path.exists(newpath):
            os.makedirs(newpath)

        import time
        t1 = time.time()
        # NOTE: To speed up, use gp.prerank instead with your own ranked list.
        # https://gseapy.readthedocs.io/en/latest/run.html
        res = gp.gsea(data=bdata.to_df().T, # row -> genes, column-> samples
            gene_sets=set,
            cls=bdata.obs.stim,
            permutation_num=1000,
            permutation_type='phenotype',
            outdir=setlabel+'/'+label,
            method='s2n', # signal_to_noise
            threads= 2,
            seed=7) # remember to set seed to get fixed seed
        t2=time.time()
        print(t2-t1)

        #res.res2d.to_csv(setlabel+'/'+label+'/gsea_seed7.out.txt', sep="\t", index=False)
        res.res2d.head(10)
        res.ranking.shape
        ## Heatmap of gene expression
        i = 1
        genes = res.res2d.Lead_genes.iloc[i].split(";")
        ax = gp.heatmap(df = res.heatmat.loc[genes],
                   z_score=None,
                   title=res.res2d.Term.iloc[i],
                   figsize=(6,5),
                   cmap=plt.cm.viridis,
                   xticklabels=False,ofname=setlabel+'/'+label+'/heatmap.png')

        ## GSEA plot
        term = res.res2d.Term
        # gp.gseaplot(res.ranking, term=term[i], **res.results[term[i]])
        axs = res.plot(terms=term[:5],ofname=setlabel+'/'+label+'/gseaplot.png')
        axs = res.plot(terms=term[:10],ofname=setlabel+'/'+label+'/gseaplot2.png')
        axs = res.plot(terms=term[:20],ofname=setlabel+'/'+label+'/gseaplot3.png')

        ## Dotplot
        # https://gseapy.readthedocs.io/en/master/gseapy_example.html#GSEA-Example
        from gseapy import dotplot
        # to save your figure, make sure that ``ofname`` is not None
        ax = dotplot(res.res2d,
                     column="FDR q-val",
                     title=set,
                     cmap=plt.cm.viridis,
                     size=5,
                     figsize=(4,5), cutoff=1,ofname=setlabel+'/'+label+'/gseaDotplot.png')

        # https://gseapy.readthedocs.io/en/master/singlecell_example.html
        # DEG Analysis
        sc.tl.rank_genes_groups(bdata,
                                groupby='stim',
                                use_raw=False,
                                layer='lognorm',
                                method='wilcoxon',
                                groups=[label2],
                                reference=label1)

        bdata.X.max() # already log1p
        sc.pl.rank_genes_groups(bdata, n_genes=25, sharey=False,save=setlabel+'_'+label+'.png')

        # get deg result
        result = bdata.uns['rank_genes_groups']
        groups = result['names'].dtype.names
        degs = pd.DataFrame(
            {group + '_' + key: result[key][group]
            for group in groups for key in ['names','scores', 'pvals','pvals_adj','logfoldchanges']})

        degs.head()
        degs.shape

        # Over-representation analysis using Enrichr
        # subset up or down regulated genes
        degs_sig = degs[degs[degs.columns[3]] < 0.05]
        degs_up = degs_sig[degs_sig[degs_sig.columns[4]] > 0]
        degs_dw = degs_sig[degs_sig[degs_sig.columns[4]] < 0]

        degs_up.shape
        degs_dw.shape

        # Enricr API
        enr_up = gp.enrichr(degs_up[degs_up.columns[0]],
                            gene_sets=set,
                            outdir=setlabel+'/'+label+'/Up')

        # trim (go:...)
        enr_up.res2d.Term = enr_up.res2d.Term.str.split(r" \(GO").str[0]
        # dotplot
        gp.dotplot(enr_up.res2d,
            figsize=(3,5), title="Up", cmap = plt.cm.autumn_r,
            ofname=setlabel+'/'+label+'/dotplotUp.png')
        plt.show()

        enr_up.res2d['-log10(Adjusted P-value)']=-np.log10(enr_up.res2d['Adjusted P-value'])
        enr_up.res2d=enr_up.res2d.sort_values(by=['Adjusted P-value'])
        gp.dotplot(enr_up.res2d, column='Combined Score',x='-log10(Adjusted P-value)',
            figsize=(3,5), title=label2, cmap = plt.cm.autumn_r,
            ofname=setlabel+'/'+label+'/dotplotUp2.png')
        plt.show()

        enr_dw = gp.enrichr(degs_dw[degs_dw.columns[0]],
                            gene_sets=set,
                            outdir=setlabel+'/'+label+'/Down')
        enr_dw.res2d.Term = enr_dw.res2d.Term.str.split(r" \(GO").str[0]
        gp.dotplot(enr_dw.res2d,
                   figsize=(3,5),
                   title="Down",
                   cmap = plt.cm.winter_r,
                   size=5,
                   ofname=setlabel+'/'+label+'/dotplotDown.png')
        plt.show()

        enr_dw.res2d['-log10(Adjusted P-value)']=-np.log10(enr_dw.res2d['Adjusted P-value'])
        enr_dw.res2d=enr_dw.res2d.sort_values(by=['Adjusted P-value'])
        gp.dotplot(enr_dw.res2d, column='Combined Score',x='-log10(Adjusted P-value)',
            figsize=(3,5), title=label1, cmap = plt.cm.autumn_r,
            ofname=setlabel+'/'+label+'/dotplotDown2.png')
        plt.show()


        # concat results
        enr_up.res2d['UP_DW'] = "UP"
        enr_dw.res2d['UP_DW'] = "DOWN"
        enr_up.res2d.sort_values(by=['Adjusted P-value'])[enr_up.res2d.columns[3:9]]
        enr_res = pd.concat([enr_up.res2d.sort_values(by=['Adjusted P-value']).head(), enr_dw.res2d.sort_values(by=['Adjusted P-value']).head()])

        from gseapy.scipalette import SciPalette
        sci = SciPalette()
        NbDr = sci.create_colormap()
        # NbDr
        # display multi-datasets
        ax = gp.dotplot(enr_res,figsize=(3,5),
                        x='UP_DW',
                        x_order = ["UP","DOWN"],
                        title=setlabel,
                        cmap = NbDr.reversed(),
                        size=3,
                        show_ring=True,ofname=setlabel+'/'+label+'/dotplot.png')
        #ax.set_xlabel("")
        plt.show()
        ax = gp.dotplot(enr_res,figsize=(3,5),
                        x='UP_DW',
                        x_order = ["UP","DOWN"],
                        title=setlabel,
                        cmap = NbDr.reversed(),
                        size=3,
                        ofname=setlabel+'/'+label+'/dotplot2.png')

        ax = gp.barplot(enr_res, figsize=(3,5),
                        group ='UP_DW',
                        title =setlabel,
                        color = ['b','r'],ofname=setlabel+'/'+label+'/barplot.png')


        # Network Visualization
        import networkx as nx
        res.res2d.head()
        nodes, edges = gp.enrichment_map(res.res2d)
        nodes.head()
        edges.head()



        # build graph
        G = nx.from_pandas_edgelist(edges,
                            source='src_idx',
                            target='targ_idx',
                            edge_attr=['jaccard_coef', 'overlap_coef', 'overlap_genes'])

        # Add missing node if there is any
        for node in nodes.index:
            if node not in G.nodes():
                G.add_node(node)


        fig, ax = plt.subplots(figsize=(8, 8))

        # init node cooridnates
        pos=nx.layout.spiral_layout(G)
        #node_size = nx.get_node_attributes()
        # draw node
        nx.draw_networkx_nodes(G,
                               pos=pos,
                               cmap=plt.cm.RdYlBu,
                               node_color=list(nodes.NES),
                               node_size=list(nodes.Hits_ratio *1000))
        # draw node label
        nx.draw_networkx_labels(G,
                                pos=pos,
                                labels=nodes.Term.to_dict())
        # draw edge
        edge_weight = nx.get_edge_attributes(G, 'jaccard_coef').values()
        nx.draw_networkx_edges(G,
                               pos=pos,
                               width=list(map(lambda x: x*10, edge_weight)),
                               edge_color='#CDDBD4')
        plt.show()
        plt.savefig(setlabel+'/'+label+'/network.png', dpi=300,bbox_inches='tight')
