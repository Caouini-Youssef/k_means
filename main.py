import random

import matplotlib
import pandas
import math
import numpy
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
matplotlib.use('macosx')



def euclidian_distance(centroid, row):
    index_df_map = {0: 'sch9/wt',
                    1: 'ras2/wt',
                    2: 'tor1/wt'}
    distance = 0.0
    for i in range(0, len(centroid)):
        tmp = row[index_df_map[i]]
        distance += (centroid[i] - row[index_df_map[i]]) ** 2
    return math.sqrt(distance)


def define_centroids(centroids, dataset):
    for index, row in dataset.iterrows():
        distances = []
        for centroid in range(0, len(centroids)):
            distances.append([euclidian_distance(centroids[centroid], row), centroid])
            distances.sort()
            row['centroid'] = distances[0][1]

def calculate_centroids(centroids, dataset):
    new_centroids = []
    grouped_dataset = dataset.groupby('centroid')
    for centroid in range(0, len(centroids)):
        group = grouped_dataset.get_group(centroid)
        mean = group.mean(axis=0)
        new_centroids.append([mean['sch9/wt'], mean['ras2/wt'], mean['tor1/wt']])
    return new_centroids

def k_means(k):
    dataset = pandas.read_excel('Longotor1delta.xls')
    dataset = dataset.drop(columns=['Public ID', 'Gene', 'Gene description'])
    previous_centroids = []
    centroids = []
    for i in range(0, k):
        previous_centroids.append([0,0,0])
        rand = random.randint(0, len(dataset.index))
        centroids.append([dataset.at[rand, 'sch9/wt'], dataset.at[rand, 'ras2/wt'], dataset.at[rand, 'tor1/wt']])
    while previous_centroids != centroids:
        define_centroids(centroids, dataset)
        previous_centroids = centroids
        centroids = calculate_centroids(centroids, dataset)
    dataset.to_csv('result.csv')

    fig = plt.figure(figsize=(10, 10))
    ax = plt.axes(projection='3d')
    ax.grid()
    grouped_dataset = dataset.groupby('centroid')
    index_c_map = {
        0: 'r',
        1: 'b',
        2: 'g',
        3: 'c',
        4: 'm',
        5: 'y'
    }
    for centroid in range(0, len(centroids)):
        group = grouped_dataset.get_group(centroid)
        ax.scatter(group['sch9/wt'], group['ras2/wt'], group['tor1/wt'], c=index_c_map[centroid], s=50)
    ax.set_title('Yeast Clustering')

    # Set axes label
    ax.set_xlabel('sch9/wt', labelpad=20)
    ax.set_ylabel('ras2/wt', labelpad=20)
    ax.set_zlabel('tor1/wt', labelpad=20)

    plt.show()

if __name__ == '__main__':
    k_means(6)
    #for i in range(1, 6):
        #k_means(i)