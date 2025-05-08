import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, CheckButtons
from scipy.signal import butter, filtfilt

# Шумова генерація
def create_noise(mu, var, time_arr):
    return np.random.normal(mu, np.sqrt(var), len(time_arr))

# Синусоїда
def make_wave(a, f, p, timeline):
    return a * np.sin(2 * np.pi * f * timeline + p)

# Сигнал із шумом
def noisy_wave(a, f, p, noise_data, timeline):
    return make_wave(a, f, p, timeline) + noise_data

# Фільтрація Баттерворта
def apply_lowpass(data, cutoff, rate):
    b, a = butter(N=4, Wn=cutoff / (0.5 * rate), btype='low')
    return filtfilt(b, a, data)

# Оновлення вмісту графіка
def refresh_plot():
    show_noise = cb_show_noise.get_status()[0]
    show_filter = cb_show_filtered.get_status()[0]

    signal_clean = make_wave(sl_a.val, sl_f.val, sl_p.val, t_values)
    line_clean.set_ydata(signal_clean)

    signal_noisy = noisy_wave(sl_a.val, sl_f.val, sl_p.val, noise_arr, t_values)
    line_noisy.set_ydata(signal_noisy)
    line_noisy.set_visible(show_noise)

    filtered_signal = apply_lowpass(signal_noisy, sl_cut.val, sl_rate.val)
    line_filtered.set_ydata(filtered_signal)
    line_filtered.set_visible(show_filter)

    fig.canvas.draw_idle()

# Оновити шум при зміні значень
def regenerate_noise(_=None):
    global noise_arr
    noise_arr = create_noise(sl_mean.val, sl_var.val, t_values)
    refresh_plot()

# Скидання слайдерів і чекбоксів
def reset_all(event):
    for sld in all_sliders:
        sld.reset()
    if cb_show_noise.get_status()[0]:
        cb_show_noise.set_active(0)
    if cb_show_filtered.get_status()[0]:
        cb_show_filtered.set_active(0)
    refresh_plot()

# Початкові параметри
a_init = 1.0
f_init = 1.0
p_init = 0.0
mean_init = 0.0
var_init = 0.1
cut_init = 2.0
rate_init = 100

t_values = np.linspace(0, 10, 1000)
noise_arr = create_noise(mean_init, var_init, t_values)

# Побудова кривих
sig_clean = make_wave(a_init, f_init, p_init, t_values)
sig_noisy = noisy_wave(a_init, f_init, p_init, noise_arr, t_values)
sig_filtered = apply_lowpass(sig_noisy, cut_init, rate_init)

# Графік
fig, ax = plt.subplots(figsize=(12, 6))
plt.subplots_adjust(left=0.2, bottom=0.52)
ax.set_ylim(-2, 2)

line_noisy, = ax.plot(t_values, sig_noisy, color='goldenrod')
line_clean, = ax.plot(t_values, sig_clean, color='steelblue')
line_filtered, = ax.plot(t_values, sig_filtered, color='crimson')

ax.set_title("Сигнал: синусоїда з шумом та фільтрацією")
ax.legend(["Шумовий сигнал", "Чистий сигнал", "Після фільтра"], loc='upper right')
line_noisy.set_visible(False)
line_filtered.set_visible(False)

# Слайдери
sl_a = Slider(plt.axes([0.2, 0.42, 0.65, 0.03]), 'Амплітуда', 0.0, 2.0, valinit=a_init)
sl_f = Slider(plt.axes([0.2, 0.37, 0.65, 0.03]), 'Частота', 0.0, 5.0, valinit=f_init)
sl_p = Slider(plt.axes([0.2, 0.32, 0.65, 0.03]), 'Фаза', 0.0, 2*np.pi, valinit=p_init)
sl_mean = Slider(plt.axes([0.2, 0.27, 0.65, 0.03]), 'Середнє шуму', -1.0, 1.0, valinit=mean_init)
sl_var = Slider(plt.axes([0.2, 0.22, 0.65, 0.03]), 'Дисперсія шуму', 0.0, 1.0, valinit=var_init)
sl_cut = Slider(plt.axes([0.2, 0.17, 0.65, 0.03]), 'Зріз частоти', 0.1, 10.0, valinit=cut_init)
sl_rate = Slider(plt.axes([0.2, 0.12, 0.65, 0.03]), 'Дискр. частота', 25.0, 200.0, valinit=rate_init)

all_sliders = [sl_a, sl_f, sl_p, sl_mean, sl_var, sl_cut, sl_rate]

# Чекбокси
cb_show_noise = CheckButtons(plt.axes([0.1, 0.04, 0.18, 0.045]), ['Шум'], [False])
cb_show_filtered = CheckButtons(plt.axes([0.42, 0.04, 0.22, 0.045]), ['Фільтр'], [False])

# Кнопка
btn_reset = Button(plt.axes([0.8, 0.04, 0.1, 0.045]), 'Скинути')
btn_reset.on_clicked(reset_all)

# Обробники подій
for sld in [sl_a, sl_f, sl_p, sl_cut, sl_rate]:
    sld.on_changed(lambda v: refresh_plot())
sl_mean.on_changed(regenerate_noise)
sl_var.on_changed(regenerate_noise)
cb_show_noise.on_clicked(lambda lbl: refresh_plot())
cb_show_filtered.on_clicked(lambda lbl: refresh_plot())

plt.show()
