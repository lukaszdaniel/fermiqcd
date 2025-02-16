import re

with open("fermiqcd.cpp", "r") as file:
    data = file.read()

regex = re.compile(r'\("(.*?)","(.*?)","?(.*?)"?\)')
d = {}
for item in regex.findall(data):
    d[item[0]] = d.get(item[0], []) + [(item[1], item[2])]

regex = re.compile(r'have\("(.*?)"\)')
for item in regex.findall(data):
    d[item] = d.get(item, [])

for key in sorted(d):
    print(f'    "{key}\\n"')
    for a, b in d[key]:
        print(f'    "   {a}={b}\\n"')
