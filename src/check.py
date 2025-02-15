name = 'check'
print(f"name = {name}")
print(f"name = {name[0:-1]}")

qwe = [1, 20000]
head = f"№   | score | polar| carbon| "
print(head)
check = "str"
for i, q in enumerate(qwe, start=1):
    numb = str(i).ljust(4)
    value = str(q).ljust(4)
    text = f"{numb}| {value}|{check}"
    print(text)

q = 0.234

q = str[q]
print(q, )
