import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from scipy.stats import wilcoxon
import numpy as np
import sys

# supprime les avertissements de copies de dataframe
# TODO : comprendre et supprimer correctement cet avertissement
pd.options.mode.chained_assignment = None
pd.set_option('display.max_columns', None)

IN = "results/iadhore/"
OUT = "results/python/"
if len(sys.argv) >= 2 and sys.argv[1] == '1' :
    IN = "results_test/iadhore/"
    OUT = "results_test/python/"

WINDOW_SIZE = 90 # nombre de gènes dans la fenêtre glissante
ANCHORS_MIN = 700 # nombre de gènes similaires minimum entre deux chromosomes homologues
ALPHA = 0.05 # risque alpha=5% pour le test statistique


"""
Lecture et traitement de la table des multiplicons :
id	genome_x	list_x	parent	genome_y	list_y	level	number_of_anchorpoints	profile_length	begin_x	end_x	begin_y	end_y	is_redundant
"""
df_multiplicons = pd.read_csv(IN + "multiplicons.txt", sep='\t', index_col='id')
#print(df_multiplicons[:3], "\nlen:", len(df_multiplicons))

"""Affichage du tableau du nbr de gènes homologues entre chaque chromosome. """
tmp = pd.concat( [ df_multiplicons['list_x'], df_multiplicons['list_y'] ] )
chromosomes = list( pd.unique(tmp)) 
chromosomes = sorted(chromosomes)
#print(chromosomes, "\n")

table_nb_anchors = [[0 for x in range(len(chromosomes))] for x in range(len(chromosomes))]  # table de 0 de taille n*n 
#print(table_nb_anchors)

# remplissage du nombre de points d'ancres entre chaque chromosome
for i in range(len(chromosomes)) :
    for j in range(i, len(chromosomes)) :
        chr1 = chromosomes[i]
        chr2 = chromosomes[j]
        
        rows = df_multiplicons[ ((df_multiplicons.list_x == chr1) & (df_multiplicons.list_y == chr2)) |
            ((df_multiplicons.list_x == chr2) & (df_multiplicons.list_y == chr1)) ]
        #print(rows[['list_x', 'list_y', 'number_of_anchorpoints']])
        anchors = rows['number_of_anchorpoints']
        total = anchors.values.sum()

        table_nb_anchors[i][j] = total
        table_nb_anchors[j][i] = total

#print((table_nb_anchors))

"""Affichage de la heat map du nombre de gènes semblables par paire de chromosomes"""
heat = go.Heatmap(z=table_nb_anchors,
                x=chromosomes,
                y=chromosomes,
                xgap=2, ygap=2
                )
layout = go.Layout(title_text="Heat map of anchor points in chromosomes", title_x=0.5, 
                width=600, height=600,
                xaxis_showgrid=False,
                yaxis_showgrid=False,
                yaxis_autorange='reversed'
                )
   
fig=go.Figure(data=[heat], layout=layout)   
fig.update_xaxes(side="top")

fig.write_html(OUT + "anchors_chr_pair.html")
#fig.show()



"""Creation de la df des triplets de chromosomes dont le nbr d'ancres est superieur à un seuil"""
MDchr = chromosomes[:17]
PPchr = chromosomes[17:]
#print(MDchr, "\n", PPchr, sep='')

df_triplets = pd.DataFrame(columns=['PP', 'MD1', 'MD2', 'anchorpoints_1', 'anchorpoints_2'])


"""Remplissage de la df grâce à la table du nombre d'ancres (gènes dupliqués)"""
def remplissage_df_triplets(table_nb_anchors, chromosomes, df) :

    for i in range(len(chromosomes)) : # parcourt triangulaire du tableau des ancres entre chromosomes
        for j in range(i, len(chromosomes)) :

            # on conserve uniquement les matchs entre 2 espèces différentes, et supérieur au seuil 
            if table_nb_anchors[i][j] >= ANCHORS_MIN and chromosomes[i][:2] != chromosomes[j][:2] :

                PP = chromosomes[i] if chromosomes[i][:2] == "Pp" else chromosomes[j]
                MD1 = chromosomes[i] if chromosomes[i][:2] == "Md" else chromosomes[j]
                # on utilise ces indices, sinon on ne sait pas si PP ou MD est en i ou j
                iMD1 = i if chromosomes[i][:2] == "Md" else j # l'indice du chromosome MD1 soit i soit j
                iPP = i + j - iMD1 # l'indice du chromosome PP soit j soit i

                # on trouve le 3e chromosome du triplet
                for k in range(len(chromosomes)) :
                    # 3e chromosome différent des 2 autres, et ses 2 matchs ont plus que le nbr minimum d'ancres, et il est de la même espèce que MD1
                    if (k != i and k != j and table_nb_anchors[i][k] >= ANCHORS_MIN and table_nb_anchors[j][k] >= ANCHORS_MIN and chromosomes[iMD1][:2] == chromosomes[k][:2]) :
                        MD2 = chromosomes[k]
                        MD = [MD1, MD2]
                        MD.sort() # trie MD1 et MD2 par ordre alphabetique, pour simplifier la suppression des triplets en doublons après

                        # crée la df d'une ligne et l'ajoute à df_triplets
                        row = pd.DataFrame({ 'PP':[PP], 
                                             'MD1':[MD[0]], 
                                             'MD2':[MD[1]], 
                                             'anchorpoints_1':[table_nb_anchors[iPP][iMD1]], 
                                             'anchorpoints_2':[table_nb_anchors[iPP][k]] })
                        df = pd.concat([df, row], ignore_index=True, axis=0)

    df.drop_duplicates(subset=['PP', 'MD1', 'MD2'], inplace=True, ignore_index=True) # suppression des triplets en double
    df.sort_values('PP', inplace=True, ignore_index=True) # tri par valeurs de chromosome PP
    return df


df_triplets = remplissage_df_triplets(table_nb_anchors, chromosomes, df_triplets)
print("chromosomes triplets :\n", df_triplets, "\n")


"""
Lecture et traitement de la table des multiplicon_pairs :
id	multiplicon		gene_x	gene_y	chr_x chr_y
"""
df_pairs = pd.read_csv(IN + "multiplicon_pairs_modified.txt", sep='\t', index_col='id')
df_pairs.drop(['code'], axis=1, inplace=True) # supprime la colonne 'code'

# jointure pour avoir les noms de chromosome de chaque gene
df_pairs_chr = df_pairs.join(df_multiplicons[['list_x', 'list_y']], on='multiplicon')
df_pairs_chr.rename(columns={'list_x':'chr_x', 'list_y':'chr_y'}, inplace=True)
#print(len(df_pairs_chr))

"""
Construction de la liste des triplets de genes chez les chromosomes du triplet
gene_PP          nb_MD1 nb_MD2 norm_MD1 norm_MD2
Prupe.3G011400.1  1      1      1        1
Prupe.3G011500.1  0      1      0        1
Prupe.3G011600.1  2      1      1        0.5
Prupe.3G011700.1  1      4      0.25     1
"""
def make_df_genes_triplet(PP, MD1, MD2, df_pairs_chr) :
    
    # récupération les paires de gènes impliquant au moins un gène de PP
    df_PPx = df_pairs_chr[ (df_pairs_chr.chr_x == PP) ] # le gène PP est dans le colonne chr_x
    df_PPy = df_pairs_chr[ (df_pairs_chr.chr_y == PP) ] # le gène PP est dans le colonne chr_y
    # échange de nom 2 colonnes d'une des deux df pour que tous les gènes du chr PP soient dans la meme colonne y
    df_PPx.rename(columns={'gene_x':'gene_y', 'gene_y':'gene_x', 'chr_x':'chr_y', 'chr_y':'chr_x'}, inplace=True)
    # concaténation des df l'une au-dessus de l'autre car elles ont les mêmes colonnes 
    df_tmp = pd.concat([df_PPy, df_PPx])

    # sélection des paires dont le MD correspond à l'un des 2 chr donnés en arguments 
    df_MD1 = df_tmp[ (df_tmp.chr_x == MD1) ]
    df_MD2 = df_tmp[ (df_tmp.chr_x == MD2) ]
    
    # regrouper par gene PP et compter pour chaque groupe le nbr de genes MD par gene PP
    # au sein d'une df, on n'a que des genes PP uniques
    df_MD1 = df_MD1.groupby(['gene_y'])['gene_x'].count().reset_index()
    df_MD2 = df_MD2.groupby(['gene_y'])['gene_x'].count().reset_index()
    
    df_MD1.rename({'gene_y':'gene_PP', 'gene_x':'nb_MD1'}, axis=1, inplace=True)
    df_MD2.rename({'gene_y':'gene_PP', 'gene_x':'nb_MD2'}, axis=1, inplace=True)
    #print(df_MD1.columns)
    #print(df_MD1[:10])
    
    #df_MD1.to_csv(OUT + "df_MD1.csv")
    #df_MD2.to_csv(OUT + "df_MD2.csv")

    ''' OLD VERSION
    # merge des 2 df, c'est-a-dire réalise une jointure inclusive
    df_triplet = pd.merge(df_MD1, df_MD2, on=['gene_y', 'chr_y'], how='outer', indicator=True)
    df_triplet = df_triplet.filter({'gene_y', 'gene_x_x', 'gene_x_y'}) # conservation des colonnes des gènes PP MD1 et MD2
    df_triplet.rename(columns={'gene_y':'gene_PP', 'gene_x_x':'gene_MD1', 'gene_x_y':'gene_MD2'}, inplace=True) # renomme les colonnes
    df_triplet = df_triplet.reindex(columns=['gene_PP', 'gene_MD1', 'gene_MD2']) # réordonne les colonnes
    df_triplet.sort_values('gene_PP', inplace=True, ignore_index=True) # trie les lignes par PP croissant, ordonné comme sur le chromosome
    '''
    # merge les deux df sur la colonne PP, en remplissant les manquantes par des 0
    df_triplet = pd.merge(df_MD1, df_MD2, on='gene_PP', how='outer').fillna(0)
    df_triplet = df_triplet.astype({'nb_MD1':'int', 'nb_MD2':'int'}) # convertit en entiers sans erreur avec NaN
    df_triplet.sort_values('gene_PP', inplace=True, ignore_index=True) # trie les lignes par PP croissant, ordonnés comme sur le chromosome

    df_triplet = normaliser_gene_PP(df_triplet)

    #df_triplet.to_csv(OUT + "tmp_" + PP + MD1 + MD2 + ".csv", index=False)
    return df_triplet

"""
Normalise entre 0 et 1 le nb de gènes MD par PP. 
Pour cela ajout des colonnes norm_MD1 et norm_MD2 qui valent 0 si MDx vaut 0,
et valent = MDx / max( MD1, MD2 ) sinon. 
"""
def normaliser_gene_PP(df_triplet) :
    df_triplet['norm_MD1'] = df_triplet.apply(lambda x: 0 if x.nb_MD1 == 0 else x.nb_MD1 / max(x.nb_MD1, x.nb_MD2), axis=1)
    df_triplet['norm_MD2'] = df_triplet.apply(lambda x: 0 if x.nb_MD2 == 0 else x.nb_MD2 / max(x.nb_MD1, x.nb_MD2), axis=1)
    return df_triplet


"""
Création d'une df mesurant le nombre de gènes conservés à chaque itération de la fenêtre le long du chromosome :
       sum_MD1  sum_MD2   rate_MD1   rate_MD2  iteration
86         NaN      NaN        NaN        NaN         -2
87         NaN      NaN        NaN        NaN         -1
88         NaN      NaN        NaN        NaN          0
89   72.416667     67.0  80.462963  74.444444          1
90   73.416667     67.0  81.574074  74.444444          2
91   74.416667     67.0  82.685185  74.444444          3
92   73.416667     67.0  81.574074  74.444444          4
93   74.416667     67.0  82.685185  74.444444          5
Avec sum_MDx la somme du nbr normalisé de gènes MDx dans la fenêtre. 
Avec rate MDx le pourcentage du nbr de MDx, calculé en faisant sum_MDx * 100 / taille fenêtre. 
"""
def make_df_window(df_genes) :
    # crée une new df pour les comptes de la fenetre glissante
    df_window = pd.DataFrame(dtype=float)
    df_window['gene_PP'] = df_genes.gene_PP
    df_window['sum_MD1'] = df_genes.norm_MD1.rolling(WINDOW_SIZE, min_periods=WINDOW_SIZE).sum()
    df_window['sum_MD2'] = df_genes.norm_MD2.rolling(WINDOW_SIZE, min_periods=WINDOW_SIZE).sum()

    # calcule le pourcentage de sum_MD dans la colonne rate_MD
    df_window['rate_MD1'] = df_window['sum_MD1'] * 100 / WINDOW_SIZE
    df_window['rate_MD2'] = df_window['sum_MD2'] * 100 / WINDOW_SIZE

    df_window.reset_index(drop=True) # réinitialise un index commencant à 0
    df_window['iteration'] = df_window.index - WINDOW_SIZE + 2
    #df_window.set_index('iteration', inplace=True)

    print(df_window[WINDOW_SIZE-5:WINDOW_SIZE+30])
    print("len : ", len(df_window))
    return df_window


"""Affiche le graphique du taux de conservation de genes au sein de 2 chromosomes dupliqués"""
def display_graph_fractionation(df_display, triplet, test_res) :
    MD1 = triplet.get("MD1")
    MD2 = triplet.get("MD2")
    PP = triplet.get("PP")

    # renomme les colonnes pour avoir une légende compréhensible
    df = df_display.rename(columns={'rate_MD1': MD1, 'rate_MD2': MD2})
    #print(df)

    # diagramme à ligne brisée du pourcentage de conservation des gènes le long de 2 chromosomes
    fig = px.line(df, x='iteration', y=[MD1, MD2])
    fig.update_layout(xaxis_title="window (size = " + str(WINDOW_SIZE) + ") iteration along " + PP,
                    yaxis_title="genome conservation rate (%)",
                    title="Fractionation biais between " + MD1 + " and " + MD2,
                    xaxis_range=[0, len(df)],
                    yaxis_range=[-1, 101] )

    fig.add_annotation(text=str(test_res),
                    xanchor='left',
                    yanchor='bottom',
                    font={'size':17, 'color':'black'},
                    x=0, y=0, showarrow=False)

    fig.write_html(OUT + PP + "_" + MD1 + "_" + MD2 + ".html")
    #fig.show()


"""Parcourt la liste des triplets de chromosomes pour en faire des graphes de biais de fractionnement"""
def analysis_each_triplet(df_triplets) :
    for triplet in df_triplets.to_dict('records') :
       analysis_one_triplet(triplet)
       print()


"""Analyse les données de biais de fractionnement d'un triplet et réalise son test statistique"""
def analysis_one_triplet(triplet) :
    # traitement des données
    df_genes_triplet = pd.DataFrame()
    df_genes_triplet = make_df_genes_triplet(triplet.get('PP'), triplet.get('MD1'), triplet.get('MD2'), df_pairs_chr)
    df_window = make_df_window(df_genes_triplet)

    # suppression des NaN des données de pourcentage de conservation des gènes
    df_display = df_window.dropna(subset=['rate_MD1', 'rate_MD2'])

    # affichage des données
    print("\n")
    [print(key,':',value) for key, value in triplet.items()]
    #print(df_display[:20])

    ## DEBUGGAGE
    #df_display.to_csv(OUT + "df_display.csv", index=False)
    #df_tmp = df_display[df_display.duplicated(subset='gene_PP', keep=False)]
    #df_tmp.to_csv(OUT + "df_duplicated.csv", index=False)


    test_res = interpretation_test(df_display)
    display_graph_fractionation(df_display, triplet, test_res)
    print("\n")


"""Réalise et interprete le test statistique wilcoxons"""
def interpretation_test(df_display) :
    res = wilcoxon(df_display['rate_MD1'], df_display['rate_MD2'])
    print(res)
    if (res.pvalue < ALPHA) : print("TEST SIGNIFICATIF : il existe un biais de fractionnement au risque alpha=", ALPHA, ". ", sep='')
    else : print("TEST NON SIGNIFICATIF : il n'existe pas de biais de fractionnement au risque alpha=", ALPHA, ". ", sep='')
    return res


"""test pour le premier triplet de la df contenant PP comme chromosome"""
def test(df_triplets, PP) :
    triplet = (df_triplets[ df_triplets.PP == PP ]).to_dict('records')[0]
    analysis_one_triplet(triplet)


"""MAIN"""
if __name__=="__main__" :
    print("Python main: running...")
    
    # lance l'analyse de tous les triplets trouvés selon le nbr de gènes similaires
    analysis_each_triplet(df_triplets)

    # lance l'analyse d'un seul triplet pour tester
    #test(df_triplets, "Pp03")

    print("Python main: done.")