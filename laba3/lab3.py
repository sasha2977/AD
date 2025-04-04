import streamlit as st
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

def load_data(path):
    if not os.path.exists(path):
        st.error(
            "File not found! Please check the following:\n"
            "1. Ensure you run script in the root directory of the repository\n"
            "2. Ensure you have run the cells in lab2.ipynb or executed lab2.py"
        )
        return None
    else:
        data = pd.read_csv(path)
        return data

def reset_filters():і
    # значення за замовчуванням
    st.session_state["option"] = 'VCI'
    st.session_state["region"] = regions[0]
    st.session_state["year_range"] = (int(data["Year"].min()), int(data["Year"].max()))
    st.session_state["week_range"] = (int(data["Week"].min()), int(data["Week"].max()))
    st.session_state["sort_asc"] = False
    st.session_state["sort_desc"] = False

path = os.path.join('lab2', 'VHI_data', 'NOAA_Ukraine_Updated.csv')
data = load_data(path)

# float to int для уникнення помилок
data["Week"] = data["Week"].astype(int) 
data["Year"] = data["Year"].astype(int)

# словник для заміни індексів(номерів) на відповідні імена регіонів
region_names = {
    1: "Вінницька",
    2: "Волинська",
    3: "Дніпропетровська",
    4: "Донецька",
    5: "Житомирська",
    6: "Закарпатська",
    7: "Запорізька",
    8: "Івано-Франківська",
    9: "Київ",
    10: "Київcька",
    11: "Кіровоградська",
    12: "Луганська",
    13: "Львівська",
    14: "Миколаївська",
    15: "Одеська",
    16: "Полтавська",
    17: "Рівненська",
    18: "Севастополь",
    19: "Сумська",
    20: "Тернопільська",
    21: "Харківська",
    22: "Херсонська",
    23: "Хмельницька",
    24: "Черкаська",
    25: "Чернівецька",
    26: "Чернігівська",
    27: "Крим"
}

if data is not None:
    st.title('Lab2')
    col1, col2 = st.columns([2,3])
    with col1:
        regions = list(region_names.values())
        st.session_state.setdefault("option", 'VCI')
        st.session_state.setdefault("region", regions[0])
        st.session_state.setdefault("year_range", (int(data["Year"].min()), int(data["Year"].max())))
        st.session_state.setdefault("week_range", (int(data["Week"].min()), int(data["Week"].max())))
        st.session_state.setdefault("sort_asc", True)
        st.session_state.setdefault("sort_desc", False)

        # 1. Dropdown для вибору часових рядів
        options = ['VCI', 'TCI', 'VHI']
        # key - ідентифікатор віджета (selectbox) у session_state
        # Він використовується для збереження стану віджета. Автоматично зберігається в session_state
        option = st.selectbox('Select:', options, key="option")

        # 2. Dropdown для вибору області
        if st.session_state["region"] not in regions:
            st.session_state["region"] = regions[0]

        index = regions.index(st.session_state["region"])

        region = st.selectbox('Select Region:', regions, index=index, key="region")

        # 3. Slider для вибору інтервалу тижнів
        min_week = data['Week'].min()
        max_week = data['Week'].max()
        
        week_range = st.slider('Select week range:', min_week, max_week, (min_week, max_week), key="week_range")

        # 4. Slider для вибору інтервалу років
        min_year = data['Year'].min()
        max_year = data['Year'].max()

        year_range = st.slider('Select year range:', min_year, max_year, (min_year, max_year), key="year_range")

        # Фільтрація даних за обраними параметрами
        # знаходить для регіону відповідний номер, наприклад 1 для 'Вінницька'
        region_index = next((key for key, value in region_names.items() if value == st.session_state["region"]), None)
        filtered_data = data[
            (data['Year'] >= year_range[0]) &
            (data['Year'] <= year_range[1]) &
            (data['Week'] >= week_range[0]) &
            (data['Week'] <= week_range[1]) &
            (data['Region'] == region_index)
        ]

        # 8. Два checkbox для сортування даних
        sort_asc = st.checkbox('Sort data in Ascending Order', key="sort_asc")
        sort_desc = st.checkbox('Sort data in Descending Order', key="sort_desc")

        if sort_asc and sort_desc:
            st.warning("Both sort options selected. Only ascending will be applied.")
            sort_desc = False
        elif sort_asc:
            filtered_data = filtered_data.sort_values(by=st.session_state["option"], ascending=True)
        elif sort_desc:
            filtered_data = filtered_data.sort_values(by=st.session_state["option"], ascending=False)

        # 5. Button для скидання фільтрів
        st.button('Reset Filters', on_click=reset_filters)

    # Колонка 2
    with col2:
        # Tabs для таблиці та графіка з відфільтрованими даними, графіка порівнянь даних по областях
        tab1, tab2, tab3 = st.tabs(["Filtered Data Table", "Time Series Plot", "Region Comparison"])

        #  Таблиця з відфільтрованими даними
        with tab1:
            st.write(filtered_data)

        # Графік з відфільтрованими даними
        with tab2:
            plt.figure(figsize=(6, 3)) 
            sns.lineplot(x=filtered_data["Year"], y=filtered_data[option])
            plt.title(f"Time Series {option} for {region} region")
            years = sorted(filtered_data["Year"].unique())
            step = 2
            selected_years = years[::step] 
            plt.xticks(selected_years, rotation=90)
            st.pyplot(plt)

        # Графік порівняння даних по областях
        with tab3:
            comparison_data = data[
                (data['Year'] >= year_range[0]) & (data['Year'] <= year_range[1]) &
                (data['Week'] >= week_range[0]) & (data['Week'] <= week_range[1])
            ]
            # групуємо дані по регіонах, обчислюємо середнє {option}
            comparison_data_grouped = comparison_data.groupby('Region')[option].mean()
            comparison_data_grouped = comparison_data_grouped.sort_values(ascending=True)

            # Заміна індексів на відповідні назви регіонів
            comparison_data_grouped.index = comparison_data_grouped.index.map(region_names)
    
            plt.figure(figsize=(6, 3)) 
            sns.barplot(x=comparison_data_grouped.index, y=comparison_data_grouped.values, palette="coolwarm")
            plt.xticks(rotation=90)
            plt.xlabel('Region')
            plt.ylabel(f'Average {option}')
            plt.title(f"Average {option} by Region") 
            
            # Виведення графіку на сторінці
            st.pyplot(plt)