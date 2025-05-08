import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, CheckButtons

# Генерація шумового сигналу
def create_random_noise(mean_val, variance_val, timeline):
    return np.random.normal(mean_val, np.sqrt(variance_val), len(timeline))

# Побудова синусоїдальної функції
def sinusoid(a, f, ph, timeline):
    return a * np.sin(2 * np.pi * f * timeline + ph)

# Додавання шуму до синусоїди
def sinusoid_with_randomness(a, f, ph, noise_array, timeline):
    return sinusoid(a, f, ph, timeline) + noise_array

# Обробка оновлення графіка
def redraw_graph():
    is_noise_visible = noise_toggle.get_status()[0]

    curve_clean = sinusoid(sl_amplitude.val, sl_frequency.val, sl_phase.val, time_line)
    curve_no_noise.set_ydata(curve_clean)

    curve_noised = sinusoid_with_randomness(sl_amplitude.val, sl_frequency.val, sl_phase.val, noise_data, time_line)
    curve_with_noise.set_ydata(curve_noised)
    curve_with_noise.set_visible(is_noise_visible)

    fig.canvas.draw_idle()

# Оновити шум
def refresh_noise(_=None):
    global noise_data
    noise_data = create_random_noise(sl_mean.val, sl_var.val, time_line)
    redraw_graph()

# Скидання параметрів
def reset_all(event):
    for sld in all_sliders:
        sld.reset()
    if noise_toggle.get_status()[0]:
        noise_toggle.set_active(0)
    redraw_graph()

# Початкові значення
a0, f0, ph0 = 1.0, 1.0, 0.0
mean0, var0 = 0.0, 0.1
time_line = np.linspace(0, 10, 1000)

# Шумовий масив
noise_data = create_random_noise(mean0, var0, time_line)

# Побудова сигналів
initial_clean = sinusoid(a0, f0, ph0, time_line)
initial_noised = sinusoid_with_randomness(a0, f0, ph0, noise_data, time_line)

# Налаштування графіка
fig, ax = plt.subplots(figsize=(12, 6))
plt.subplots_adjust(left=0.25, bottom=0.45)
ax.set_ylim(-2, 2)
curve_with_noise, = ax.plot(time_line, initial_noised, color='tomato')
curve_no_noise, = ax.plot(time_line, initial_clean, color='navy')
ax.set_title("Синусоїда із шумовою складовою")
ax.legend(["З шумом", "Чистий сигнал"], loc='upper right')
curve_with_noise.set_visible(False)

# Слайдери
sl_amplitude = Slider(plt.axes([0.25, 0.35, 0.65, 0.03]), 'Амплітуда', 0.0, 2.0, valinit=a0)
sl_frequency = Slider(plt.axes([0.25, 0.30, 0.65, 0.03]), 'Частота', 0.0, 5.0, valinit=f0)
sl_phase = Slider(plt.axes([0.25, 0.25, 0.65, 0.03]), 'Фаза', 0.0, 2*np.pi, valinit=ph0)
sl_mean = Slider(plt.axes([0.25, 0.20, 0.65, 0.03]), 'Середнє Шуму', -1.0, 1.0, valinit=mean0)
sl_var = Slider(plt.axes([0.25, 0.15, 0.65, 0.03]), 'Дисперсія Шуму', 0.0, 1.0, valinit=var0)

all_sliders = [sl_amplitude, sl_frequency, sl_phase, sl_mean, sl_var]

# Чекбокс та кнопка
noise_toggle = CheckButtons(plt.axes([0.25, 0.05, 0.15, 0.04]), ['Шум'], [False])
reset_btn = Button(plt.axes([0.8, 0.05, 0.1, 0.04]), 'Скидання')
reset_btn.on_clicked(reset_all)

# Зв’язування подій
for sld in [sl_amplitude, sl_frequency, sl_phase]:
    sld.on_changed(lambda val: redraw_graph())
sl_mean.on_changed(refresh_noise)
sl_var.on_changed(refresh_noise)
noise_toggle.on_clicked(lambda lbl: redraw_graph())

plt.show()
