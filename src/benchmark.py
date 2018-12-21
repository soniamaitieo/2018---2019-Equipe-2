
import matplotlib.pyplot as plt


def parse_benchlist(bench_list):
    " Parse benmark.list and create dictionnary: Query: Familly"
    dict_bench = {}
    with open(bench_list ,'r') as f:
        lines = f.readlines()
        for line in lines:
            if line[0] == '-':
                line = line.rstrip('\n')
                splitted_line = line .split(': ')            
                query_name = splitted_line[0][2:]
                query_fam = splitted_line[1].split(' ')[0]
                dict_bench[query_name] = query_fam
    return( dict_bench)

        

def find_rank(path_to_foldrec , query_name):
    path_query = path_to_foldrec + query_name + ".foldrec"
    #path_query = "/home/sdv/m2bi/stieo/2018---2019-Equipe-2/results/Cohesin.foldrec"
    mydico = parse_benchlist(bench_list)
    with open(path_query ,'r') as myFile:
        lines = myFile.readlines()
    for num , line in enumerate(lines, 0):
        if mydico[query_name] in line:
            mynum = num
            break
    return((query_name , int (lines[mynum].split(' ')[0])))


L = 405*[0]
for t in range(0 , len(test)):
    el = test[t][1]
    for l in range(el, len(L)):
        L[l] += 1
        print(L)




def enrich_list(liste_tuple):
    "liste de tuple : liste_tuple = [('Cohesin', 26) , ('Machin', 2) , ('Truc', 10)]"
    L = 405*[0]
    for t in range(0 , len(liste_tuple)):
        el = liste_tuple[t][1]
        for l in range(el, len(L)):
            L[l] += 1
    return(L)


#-----------------------------------------------

query_name = "Cohesin"
find_rank(query_name)
 

#------------------------------------------------ PLOT 

dotprod_list = [10, 12 , 13, 20, 24, 32, 50, 100, 101, 105, 300, 500]
pearson_list = [10, 50, 100, 101, 105, 300, 400, 450, 475, 500 , 501, 502]


def plot_enrichmentcurve(dotprod_list, pearson_list):
    plt.grid(True)
    plt.plot(dotprod_list, "b", linewidth=0.8, marker="+", label="Dot Product")
    plt.plot(pearson_list, "g", linewidth=0.8, marker="+", label="Pearson Correlation")
    plt.axis([1, len(dotprod_list) + 1 , 0, 600])
    plt.xlabel('Rank')
    plt.ylabel('Enrichment')
    plt.legend()    
    plt.show()




#------------------ MAIN

path_to_foldrec_pearson = "/home/sdv/m2bi/stieo/2018---2019-partage/Data/outputs_ORION/"
path_to_foldrec_dotprod = "/home/sdv/m2bi/stieo/2018---2019-Equipe-2/results/"
bench_list = "/home/sdv/m2bi/stieo/2018---2019-partage/Data/Benchmark.list"


liste_tuple_dotprod = []
liste_tuple_pearson = []
mydico = parse_benchlist(bench_list)

for k in mydico.keys() :
    liste_tuple_dotprod.append(find_rank(path_to_foldrec_dotprod  , k))
    liste_tuple_pearson.append(find_rank(path_to_foldrec_pearson , k))
dotprod_list =  enrich_list(liste_tuple_dotprod) 
pearson_list =  enrich_list(liste_tuple_pearson) 
plot_enrichmentcurve(dotprod_list, pearson_list)