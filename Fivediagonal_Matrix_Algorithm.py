import numpy as np
import random
import matplotlib.pyplot as plt
import time

#Генерация матрицы

def Generate_Matrix(size):

    err=np.full(size,0.1)
    c = np.random.randint(size+1,size+20,(size,))
    d = np.random.random_sample((size,))+err
    e = np.random.random_sample((size,))+err
    b = np.random.random_sample((size,))+err
    a = np.random.random_sample((size,))+err
    f = np.random.randint(size+1,size+20,(size,))

    return (size-1,c,d,e,b,a,f)

def Algorithm(size):

    N,c,d,e,b,a,f = Generate_Matrix(size)

    start_time = time.time()

    alpha=np.zeros(N+1)
    gamma=np.zeros(N+2)
    betta=np.zeros(N)
    delta=np.zeros(N+1)

    #Начальные условия

    alpha[1]=d[0]/c[0]
    gamma[1]=f[0]/c[0]
    betta[1]=e[0]/c[0]
    delta[1]=c[1]-b[1]*alpha[1]

    alpha[2]=(d[1]-betta[1]*b[1])/delta[1]
    gamma[2]=(f[1]+b[1]*gamma[1])/delta[1]
    betta[2]=e[1]/delta[1]

    #Рекурсивный поиск коэффицинтов

    for i in range (2,N+1):
        delta[i]=c[i]-a[i]*betta[i-1]+alpha[i]*(a[i]*alpha[i-1]-b[i])
        if (i<=N-1):
            alpha[i+1]=(d[i]+betta[i]*(a[i]*alpha[i-1]-b[i]))/delta[i]
        gamma[i+1]=(f[i]-a[i]*gamma[i-1]-gamma[i]*(a[i]*alpha[i-1]-b[i]))/delta[i]
        if (i<=N-2):
            betta[i+1]=e[i]/delta[i]

    #Обратный ход

    y=np.zeros(N+1)

    y[N]=gamma[N+1]

    y[N-1]=alpha[N]*y[N]+gamma[N]

    for i in range (N-2,-1,-1):
        y[i]=alpha[i+1]*y[i+1]-betta[i+1]*y[i+2]+gamma[i+1]

    #Результаты

    # print("--- %s seconds ---" % (time.time() - start_time))

    #Проверка и результат

    # N=N+1
    # A=np.zeros((N+1,N+1))

    # d=d[:N-1]
    # e=e[:N-2]
    # b=b[1:]
    # a=a[2:]
    # A=A+np.diag(c,0)+np.diag(-d,1)+np.diag(e,2)+np.diag(-b,-1)+np.diag(a,-2)

    # print("Невязка:")
    # print(np.matmul(A,y)-f)

    # print("Размер матрицы:")
    # print(N)

    # print("Решение:")
    # print(y)

    return(time.time() - start_time)


def Test(coef,max_N):

    sizes=[]
    times=[]

    averaging=0
    for j in range(5,max_N):
        averaging+=Algorithm(j)
        if (j % 25==0):
            sizes.append(j)
            times.append(((averaging)/20)*1000000)
            averaging=0
    
    return (sizes,times)

#Тестирование на временную сложность

testing_count=10
s_list=[]
t_list=[]

for j in range(1,testing_count):
    s, t = Test(1,10000)
    s_list.append(s)
    t_list.append(t)

s=np.array(s_list[0])
t=np.array(t_list[0])

for j in range(1,len(s_list)):
    s+=np.array(s_list[j])
    t+=np.array(t_list[j])

s=s/testing_count
t=t/testing_count

#Визуализация

plt.title("Зависимость времени от размерности матрицы")
plt.xlabel("Размерность матрицы")
plt.ylabel("Время выполнения алгоритма T (мкс)")
plt.plot(s,t)
plt.show()

#Тестирование при удвоении входных данных

testing_count=5
s=[2500,5000,10000]
t_list=[]

for l in s:
    t_mid=0
    for j in range(testing_count):
        t = Algorithm(l)
        t_mid+=t
    t_list.append(t/testing_count)

print(t_list)

