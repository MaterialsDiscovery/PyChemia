import sys

data = open(sys.argv[1]).readlines()

Sections = {}

for i in range(len(data)):
    Variable = ''
    Section = ''
    if data[i][:8] == 'Variable':
        Variable = data[i].split()[1]
        for j in range(1, 5):
            if data[i + j][:8] == 'Section ':
                Section = data[i + j][8:].split('::')[0].strip()
                break
        if Section == '':
            print('Error finding section for ', Variable)

        if Section not in Sections.keys():
            print('New Section ', Section)
            Sections[Section] = [Variable]
        else:
            Sections[Section].append(Variable)

variables = open('octopus_variables.conf', 'w')
for i in sorted(Sections.keys()):
    variables.write('\n[' + i + ']\n')
    for j in sorted(Sections[i]):
        variables.write(j + '\n')
variables.close()
