import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from scipy.stats import wilcoxon
import numpy as np
import sys

# supprime les avertissements de copies de dataframe
# TODO : comprendre et supprimer correctement cet avertissement
pd.options.mode.chained_assignment = None
pd.options.display.max_rows = 999

IN = "results/iadhore/"
OUT = "results/python/"
PATH_PP = "data/PP_lst/"

if len(sys.argv) >= 2 and sys.argv[1] == '1' :
    IN = "results_test/iadhore/"
    OUT = "results_test/python/"

WINDOW_SIZE = 80 # nombre de gènes dans la fenêtre glissante, conseillé entre 50 et 100 gènes par fenêtres
ANCHORS_MIN = 200 # nombre de gènes similaires minimum entre deux chromosomes homologues, conseillé n'importe quelle valeur entre 200 et 700
ALPHA = 0.05 # risque alpha=5% pour le test statistique

# si plus de MIN_WINDOWS gènes de PP à suivre ont un nombre de gène MD conservés inférieur ou égale à MIN_RATE, alors on exclue ces données de l'analyse statistique, car c'est un fragment entier du chromosome qui manque, ce n'est pas un effet du fractionation biais
NO_SYNTENY_MIN_WINDOWS = 200 # nombre de gènes dont la conservation est inférieur au seuil, à partir duquel on considère que le fragment n'est pas synténique, conseillé entre 100 et 500
NO_SYNTENY_MIN_RATE = 0 # borne inférieure exclue des gènes en synténie, fortement conseillé 0 (inférieur à 1)


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
gene_PP          nb_MD1 nb_MD2
Prupe.3G011400.1  1      1    
Prupe.3G011500.1  0      1    
Prupe.3G011700.1  1      4    
"""
def make_df_genes_triplet(PP, MD1, MD2, df_pairs_chr) :
    
    # récupération les paires de gènes impliquant au moins un gène de PP
    df_PPx = df_pairs_chr[ (df_pairs_chr.chr_x == PP) ] # le gène PP est dans le colonne chr_x
    df_PPy = df_pairs_chr[ (df_pairs_chr.chr_y == PP) ] # le gène PP est dans le colonne chr_y
    # échange de nom 2 colonnes d'une des deux df pour que tous les gènes du chr PP soient dans la meme colonne y
    df_PPx.rename(columns={'gene_x':'gene_y', 'gene_y':'gene_x', 'chr_x':'chr_y', 'chr_y':'chr_x'}, inplace=True)
    # concaténation des df l'une au-dessus de l'autre car elles ont les mêmes colonnes 
    df_tmp = pd.concat([df_PPy, df_PPx])
    #print(df_tmp)

    # sélection des paires dont le MD correspond à l'un des 2 chr donnés en arguments 
    df_tmp_MD1 = df_tmp[ (df_tmp.chr_x == MD1) ]
    df_tmp_MD2 = df_tmp[ (df_tmp.chr_x == MD2) ]
    
    df_tmp_MD1.rename({'gene_y':'gene_PP', 'gene_x':'gene_MD1'}, axis=1, inplace=True)
    df_tmp_MD2.rename({'gene_y':'gene_PP', 'gene_x':'gene_MD2'}, axis=1, inplace=True)

    # regroupe les gènes MD en listes pour avoir un PP unique par ligne
    df_MD1 = df_tmp_MD1.groupby('gene_PP')['gene_MD1'].apply(list).reset_index(name='gene_MD1')
    df_MD2 = df_tmp_MD2.groupby('gene_PP')['gene_MD2'].apply(list).reset_index(name='gene_MD2')

    # regrouper par gene PP et compter pour chaque groupe le nbr de genes MD par gene PP
    # au sein d'une df, on n'a que des genes PP uniques
    df_nbMD1 = df_tmp_MD1.groupby(['gene_PP'])['gene_MD1'].count().reset_index(name='nb_MD1')
    df_nbMD2 = df_tmp_MD2.groupby(['gene_PP'])['gene_MD2'].count().reset_index(name='nb_MD2')

    # merge les quatre df sur la colonne PP, en remplissant les manquantes par des 0
    df_triplet = pd.merge(df_MD1, df_MD2, on='gene_PP', how='outer').fillna(0)
    df_triplet = pd.merge(df_triplet, df_nbMD1, on='gene_PP', how='outer').fillna(0)
    df_triplet = pd.merge(df_triplet, df_nbMD2, on='gene_PP', how='outer').fillna(0)

    df_triplet = df_triplet.astype({'nb_MD1':'int', 'nb_MD2':'int'}) # convertit en entiers sans erreur avec NaN
    df_triplet.sort_values('gene_PP', inplace=True, ignore_index=True) # trie les lignes par PP croissant, ordonnés comme sur le chromosome

    return df_triplet


"""Ajout des gènes manquants de PP, avec des valeurs de 0 pour les MD correspondants"""
def add_every_PP(df_triplet, PP) :
    # liste les gènes du génome complet du chromosome du triplet
    df_PP = pd.read_csv(PATH_PP + PP + ".lst", usecols=[0], names=['gene_PP'])
    # supprime l'orientation (+ ou -) à la fin de chaque gene
    df_PP['gene_PP'] = df_PP.apply( lambda x: x.gene_PP[: len(x.gene_PP)-1 ] , axis=1)
    
    # concatène les df l'une au-dessus de l'autre puis supprime les duplicata de genes PP en gardant la premiere occurrence car c'est celle de df_triplet
    # on obtient bien le meme nbr de lignes que de genes dans le fichier .lst
    df_triplet = pd.concat([df_triplet, df_PP]).drop_duplicates(subset=['gene_PP'], keep='first').reset_index(drop=True)
    df_triplet.fillna(0, inplace=True) # remplie les NaN avec des 0
    # ordonne par gene PP croissant
    df_triplet = df_triplet.sort_values(by=['gene_PP'], ignore_index=True)

    return df_triplet

"""
Normalise entre 0 et 1 le nb de gènes MD par PP. 
Pour cela ajout des colonnes norm_MD1 et norm_MD2 qui valent 0 si MDx vaut 0,
et valent = MDx / max( MD1, MD2 ) sinon. 
gene_PP          nb_MD1 nb_MD2 norm_MD1 norm_MD2
Prupe.3G011400.1  1      1      1        1
Prupe.3G011500.1  0      1      0        1
Prupe.3G011700.1  1      4      0.25     1
"""
def normaliser_gene_PP(df_triplet) :
    df_triplet['norm_MD1'] = df_triplet.apply(lambda x: 0 if x.nb_MD1 == 0 else x.nb_MD1 / max(x.nb_MD1, x.nb_MD2), axis=1)
    df_triplet['norm_MD2'] = df_triplet.apply(lambda x: 0 if x.nb_MD2 == 0 else x.nb_MD2 / max(x.nb_MD1, x.nb_MD2), axis=1)
    return df_triplet

"""
Création d'une df mesurant le nombre de gènes conservés à chaque itération de la fenêtre le long du chromosome :
       sum_MD1  sum_MD2   rate_MD1   rate_MD2  iteration
87         NaN      NaN        NaN        NaN         -1
88         NaN      NaN        NaN        NaN          0
89   72.416667     67.0  80.462963  74.444444          1
90   73.416667     67.0  81.574074  74.444444          2
91   74.416667     67.0  82.685185  74.444444          3
Avec sum_MDx la somme du nbr normalisé de gènes MDx dans la fenêtre. 
Avec rate MDx le pourcentage du nbr de MDx, calculé en faisant sum_MDx * 100 / taille fenêtre. 
"""
def make_df_window(df_triplet) :
    # crée une new df pour les comptes de la fenetre glissante
    df_window = pd.DataFrame(dtype=float)
    df_window['gene_PP'] = df_triplet.gene_PP
    df_window['sum_MD1'] = df_triplet.norm_MD1.rolling(WINDOW_SIZE, min_periods=WINDOW_SIZE, center=True).sum()
    df_window['sum_MD2'] = df_triplet.norm_MD2.rolling(WINDOW_SIZE, min_periods=WINDOW_SIZE, center=True).sum()

    # calcule le pourcentage de sum_MD dans la colonne rate_MD
    df_window['rate_MD1'] = df_window['sum_MD1'] * 100 / WINDOW_SIZE
    df_window['rate_MD2'] = df_window['sum_MD2'] * 100 / WINDOW_SIZE

    # arrondis les pourcentage à 3 décimales
    df_window.rate_MD1 = df_window.rate_MD1.astype(float).round(decimals = 3)
    df_window.rate_MD2 = df_window.rate_MD2.astype(float).round(decimals = 3)

    df_window.reset_index(drop=True) # réinitialise un index commencant à 0
    df_window['iteration'] = df_window.index - WINDOW_SIZE // 2 + 1

    #print(df_window[WINDOW_SIZE-5:WINDOW_SIZE+30])
    #print("len : ", len(df_window))
    return df_window


"""Affiche le graphique du taux de conservation de genes au sein de 2 chromosomes dupliqués"""
def display_graph_fractionation(df_display, triplet, test_res, df_synteny) :
    MD1 = triplet.get("MD1")
    MD2 = triplet.get("MD2")
    PP = triplet.get("PP")

    # renomme les colonnes pour avoir une légende compréhensible
    df = df_display.rename(columns={'rate_MD1': MD1, 'rate_MD2': MD2})

    # diagramme à ligne brisée du pourcentage de conservation des gènes le long de 2 chromosomes
    fig = px.line(df, x=df.index, y=[MD1, MD2])
    
    fig.update_layout(xaxis_title="window (size = " + str(WINDOW_SIZE) + ") iteration along " + PP,
                    yaxis_title="genome conservation rate (%)",
                    title="Fractionation biais between " + MD1 + " and " + MD2 + " (compared to " + PP + ")",
                    xaxis_range=[0, len(df)],
                    yaxis_range=[-1, 101] )

    fig.add_annotation(text=str(test_res),
                    xanchor='left',
                    yanchor='bottom',
                    font={'size':17, 'color':'black'},
                    x=0, y=0, showarrow=False)
    
    for index, row in df_synteny.iterrows() :
        fig.add_vrect(x0=row['debut'], x1=row['fin'], 
                    annotation_text="synténie (" + str(row['debut']) + "-" + str(row['fin']) + ")", 
                    annotation_position="top left",
                    fillcolor="green", 
                    opacity=0.2, 
                    line_width=0, 
                    layer="below")

    fig.write_html(OUT + PP + "_" + MD1 + "_" + MD2 + ".html")
    #fig.show() # ne fonctionne pas en ssh ?


""" 
Trouve les limites des blocs de synténies ou plusieurs valeurs du nb de MD1 sont inférieurs à un seuil
                        >seuil  somme   groupby+count   <seuil  rolling     diff
itération norm_MD1      tmp1    tmp2    tmp3            tmp4    synteny     limit
       1        1        1      1       2               1       1           NaN
       2        0        0      1       2               1       0           1
       3        1        1      2       4               0       0           0
       4        0        0      2       4               0       0           0
       5      0.2        0      2       4               0       0           0
       6        0        0      2       4               0       0           0
       7        1        1      3       1               1       0           0
       8        1        1      4       1               1       1          -1
       9        1        1      5       1               1       1           0
"""
def make_synteny_limits(df_display) :
    # pour MD1
    df_display['tmp_MD1'] = df_display.apply(lambda x: 0 if x.norm_MD1 <= NO_SYNTENY_MIN_RATE else 1 , axis=1)
    df_display['tmp2_MD1'] = df_display.tmp_MD1.cumsum()
    df_display['tmp3_MD1'] = df_display.groupby('tmp2_MD1')['tmp2_MD1'].transform('count')
    df_display['tmp4_MD1'] = df_display.apply(lambda x: 0 if x.tmp3_MD1 > NO_SYNTENY_MIN_WINDOWS else 1, axis=1)
    # toutes les fenetres contenant au moins une valeur hors du bloc de synténie est hors synténie
    df_display['synteny_MD1'] = df_display.tmp4_MD1.rolling(WINDOW_SIZE, min_periods=1, center=True).min()

    # idem pour MD2
    df_display['tmp_MD2'] = df_display.apply(lambda x: 0 if x.norm_MD2 <= NO_SYNTENY_MIN_RATE else 1 , axis=1)
    df_display['tmp2_MD2'] = df_display.tmp_MD2.cumsum()
    df_display['tmp3_MD2'] = df_display.groupby('tmp2_MD2')['tmp2_MD2'].transform('count')
    df_display['tmp4_MD2'] = df_display.apply(lambda x: 0 if x.tmp3_MD2 > NO_SYNTENY_MIN_WINDOWS else 1, axis=1)
    df_display['synteny_MD2'] = df_display.tmp4_MD2.rolling(WINDOW_SIZE, min_periods=1, center=True).min()

    # sont en synténie seulement les blocs ou MD1 et MD2 sont en synténie
    df_display['synteny'] = df_display.apply(lambda x: 0 if x.synteny_MD1 == 0 or x.synteny_MD2 == 0 else 1, axis=1)
    df_display['limit'] = df_display.synteny.diff()

    df_display.drop(['tmp_MD1', 'tmp2_MD1', 'tmp3_MD1', 'tmp4_MD1', 'tmp_MD2', 'tmp2_MD2', 'tmp3_MD2', 'tmp4_MD2'], axis=1, inplace=True) # supprime les colonnes temporaires
    
    df_synteny = pd.DataFrame(df_display.loc[df_display['limit'] != 0, 'limit']).reindex()
    end = len(df_display)
    if (len(df_synteny) > 1) :
        df_synteny = traiter_synteny(df_synteny, end)
    else : # tout le chromosome en synténie, donc dataframe une seule ligne : debut = 1 et fin = fin du chromosome
        df_synteny = pd.DataFrame([[1, end]], columns=['debut', 'fin'])

    return df_synteny, df_display

"""
Trouve les limites des fragments de synténie en créant la df qui indique chaque début et fin de bloc de synténie :
    debut   fin
0    346   801
1    945  1069
2   2283  3583
"""
def traiter_synteny(df_synteny, end) :
    df_synteny.reset_index(inplace=True)
    if df_synteny.iloc[1]['limit'] == -1 :
        df_synteny.at[0, 'limit'] = 1 
    if df_synteny.iloc[len(df_synteny) - 1]['limit'] == 1 :
        df_synteny.loc[len(df_synteny)] = [end, -1]

    df_debut = df_synteny[df_synteny.limit == 1].reset_index(drop=True)
    df_debut.rename(columns={'iteration':'debut'}, inplace=True)
    df_debut.drop(['limit'], axis=1, inplace=True)
    df_fin = df_synteny[df_synteny.limit == -1].reset_index(drop=True)
    df_fin.rename(columns={'iteration':'fin'}, inplace=True)
    df_fin.drop(['limit'], axis=1, inplace=True)
    df_synteny = pd.concat([df_debut, df_fin], axis=1)

    return df_synteny


"""Analyse les données de biais de fractionnement d'un triplet et réalise son test statistique"""
def analysis_one_triplet(triplet) :
    PP = triplet.get('PP')
    MD1 = triplet.get('MD1')
    MD2 = triplet.get('MD2')

    # crée la df du nbr de MD par triplets de gènes
    df_genes_triplet = make_df_genes_triplet(PP, MD1, MD2, df_pairs_chr)

    # ajoute les gènes de PP manquants
    df_genes_triplet = add_every_PP(df_genes_triplet, PP)

    # normalise les nbr de MD entre 0 et 1
    df_genes_triplet = normaliser_gene_PP(df_genes_triplet)

    # crée la df des valeurs en fenetre glissante, avec un pas=1
    df_window = make_df_window(df_genes_triplet)

    # merge les df des gènes avec celle des pourcentages, en centrant la fenetre de pourcentage sur le gène
    # et exporte ces données pour être exploitées par Tanguy
    df_data_complete = pd.merge(df_genes_triplet, df_window, on='gene_PP', how='outer')
    df_data_complete = df_data_complete.explode('gene_MD1')
    df_data_complete = df_data_complete.explode('gene_MD2')
    df_data_complete.to_csv(OUT + "gene_rate_" + PP + "_" + MD1 + "_" + MD2 + ".csv", index=False)

    # suppression des NaN des données de pourcentage de conservation des gènes
    df_display = df_window.dropna(subset=['rate_MD1', 'rate_MD2'])
    df_display = df_display.merge(df_genes_triplet[['gene_PP', 'norm_MD1', 'norm_MD2']], on='gene_PP')
    df_display = df_display.set_index('iteration', drop=True)
    #df_display.to_csv(OUT + "display_" + PP + "_" + MD1 + "_" + MD2 + ".csv")

    df_synteny, df_display = make_synteny_limits(df_display)

    # affichage des données
    [print(key,':',value) for key, value in triplet.items()]
    print("\nSYNTENY : \n", df_synteny)
    test_res = interpretation_test(df_display)
    display_graph_fractionation(df_display, triplet, test_res, df_synteny)
    print("\n==========================================\n")


"""Parcourt la liste des triplets de chromosomes pour en faire des graphes de biais de fractionnement"""
def analysis_each_triplet(df_triplets) :
    for triplet in df_triplets.to_dict('records') :
       analysis_one_triplet(triplet)


"""Réalise et interprete le test statistique wilcoxons sur les blocs de synténie"""
def interpretation_test(df_display) :
    df_test = df_display[df_display.synteny == 1]
    res = wilcoxon(df_test['rate_MD1'], df_test['rate_MD2'])
    print("\n" + str(res))
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
    
    # lance l'analyse sur tous les triplets trouvés
    analysis_each_triplet(df_triplets)

    # lance l'analyse d'un seul triplet pour tester
    #test(df_triplets, "Pp03")

    print("Python main: done.")