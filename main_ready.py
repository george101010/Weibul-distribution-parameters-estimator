import tkinter as tk
import numpy as np
from scipy.stats import weibull_min , norm


# A_ , B_  - границы поиска параметра бета
# X_ - значения в выборке
# n - всего испытаний
#  r - количество элементов в выборке
# max_iter = 50 - критерий остановки
#  tol = 0.0001 - критерий остановки
def get_params_weibull(A_ , B_  , X_ , n , r , max_iter = 50 , tol = 0.0001 ):
    def equation_beta( Beta , X_ , n , r ):

        alpha_ = r / ( (n - r) * X_[-1] **Beta + np.sum(  X_ **Beta  ) )
        beta_f =  r / Beta - (n - r) * alpha_ * X_[-1] **Beta * np.log(X_[-1]) + np.sum(  np.log(X_ ) * (1 - alpha_ * X_ **Beta) )
        return beta_f

    def get_alpha( Beta  , X_ , n , r):
        alpha_ = r / ( (n - r) * X_[-1] **Beta + np.sum(  X_ **Beta  ) )
        return alpha_

    def solver( f  ,A_ , B_  ,*args) :
        if f(A_,*args) * f( B_ ,*args) >= 0:
            raise ValueError("Корень не находится на заданном интервале.")
        max_iter = 50
        tol = 0.0001
        for i in range(max_iter):
            C_ = ( B_  +  A_ ) / 2
            if abs(f(C_,*args)) < tol:
                return C_
            elif f(A_,*args) * f(C_,*args) < 0:
                B_ = C_
            else:
                A_ = C_

        raise ValueError("Метод не сошелся за заданное число итераций.")

    solved_Beta = solver(equation_beta , A_ ,B_ , X_ , n , r )
    solved_alpha = get_alpha( solved_Beta  ,X_  , n , r)

    return { 'alpha': solved_alpha  , 'beta' : solved_Beta  }

# доверительный интервал по найденным параметрам
def conf_interval_weibull( Alpha_ , Beta_ , Sigma , X_ , r , n):

    dd_alpha = - r / ( Alpha_ * Alpha_ )
    dd_beta = - r / ( Beta_ * Beta_ ) - (n-r) * Alpha_ *X_[-1]**Beta_ * np.log(X_[-1])* np.log(X_[-1]) - Alpha_ * np.sum(X_**Beta_ * np.log(X_)* np.log(X_))
    dd_alpha_beta = - (n-r) * X_[-1]**Beta_ * np.log(X_[-1])  - np.sum(X_**Beta_ * np.log(X_) )

    Det = dd_alpha * dd_beta - dd_alpha_beta * dd_alpha_beta

    Edge_alpha = - dd_beta / Det
    Edge_beta = - dd_alpha / Det

    Norm_proc = norm.ppf( 1 -  Sigma / 2  )

    delta_alpha = Norm_proc * Edge_alpha ** 0.5
    delta_beta = Norm_proc * Edge_beta ** 0.5

    Interval_alpha =  ( Alpha_ - delta_alpha, Alpha_ + delta_alpha)
    Interval_beta = ( Beta_ - delta_beta, Beta_ + delta_beta)

    return { "alpha" :  Interval_alpha,
             "beta" : Interval_beta }

def get_display_result():
    # Получение всех данных
    #samples = np.array( list( map(float , entry_sample.get().split(','))))
    samples =  np.array( list( map(float , entry_sample.get().split(','))) )
    N_sample = float( entry_n.get() ) 
    
    right_brd =  float( entry_rb.get() ) 
    left_brd =  float( entry_lb.get() ) 
    percent = int(entry_percent.get() )

    Sigma_ = 1 - percent / 100

    params = get_params_weibull( right_brd  ,left_brd , samples , N_sample , len(samples ) , max_iter = 1000 , tol = 0.0001  )
    
    confidence_intervals = conf_interval_weibull( params['alpha'] , params['beta']  , Sigma_ , samples ,  len(samples ) , N_sample)


    OUTPUT_TEXT = f"Оцененные параметры:\nAlpha: {np.round(params['alpha'] , 16)}\nBeta: {np.round(params['beta'] , 16)}\n"
    OUTPUT_TEXT += f"Доверительный интервал:\nAlpha: [{np.round(confidence_intervals['alpha'][0] , 16)} , {np.round(confidence_intervals['alpha'][1] , 16)}]\nBeta: [{np.round(confidence_intervals['beta'][0] , 16)} , {np.round(confidence_intervals['beta'][1] , 16)}]"    
    # Отображаем текст в текстовом поле
    output_text.delete(1.0, tk.END)  # Очищаем текстовое поле
    output_text.insert(tk.END,   OUTPUT_TEXT  )
 
# Создаем основное окно
main_window = tk.Tk()
main_window.title("Распределение Вейбулла. Оценка параметров.")
main_window.config(width = 1000, height = 360 )
main_window.resizable( False, False )
# подпись - введите элементы выборки
label_sample = tk.Label( text = "Введите элементы выборки через запятую:")
label_sample.place( x = 0, y = 0+4,  width = 330 , height = 30)  
# поле ввода - элементы выборки
entry_sample = tk.Entry(main_window )
entry_sample.place(x = 331 , y =  0+4,  width = 640 , height = 30)  

# подпись - введите элементы выборки
label_n = tk.Label( text = "Введите общий размер выборки:")
label_n.place( x = 0, y = 30+4,  width = 270 , height = 30)  
# поле ввода - общее количество элементов
entry_n = tk.Entry(main_window  )
entry_n.place(x =  271, y =  30+4,  width = 40, height = 30)  
# подпись - введите границы
label_n = tk.Label( text = "Введите границы поиска параметра beta для метода половинниго деления:")
label_n.place( x = 0, y = 60+8,  width = 560 , height = 30)  
 
# поле ввода - правая граница поиска
entry_rb = tk.Entry(main_window  )
entry_rb.place(x = 560 + 4 , y = 60+8 ,  width = 40 , height = 30)  
# поле ввода - левая граница поиска
entry_lb = tk.Entry(main_window )
entry_lb.place(x =   560 + 40+8, y = 60+8 ,  width = 40, height = 30 ) 

# подпись - введите параметр пси
label_percent = tk.Label( text = "Введите значение % для доверительного интервала:")
label_percent.place( x = 0, y = 90+12,  width = 380 , height = 30)   
# поле ввода - левая граница поиска
entry_percent = tk.Entry(main_window )
entry_percent.place(x =   380 + 4, y = 90+12 ,  width = 40, height = 30 ) 

 
 
# Создаем кнопку для отображения текста
button = tk.Button(main_window, text="Оценка параметров", command=  get_display_result)
button.place(x = 140 , y = 120 +16,  width = 180, height =30 )  
 
# Создаем текстовое поле для отображения информации
output_text = tk.Text(main_window )
output_text.place(x = 4 , y =  150+20,  width = 700, height = 180)  
 
# Запускаем главный цикл приложения
main_window.mainloop()
