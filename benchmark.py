#!/usr/bin/env python2
from sys import argv
from csv import writer as csv_writer
from numpy import median
from time import strftime 

def main():
    if argv[1] == '--preprocess':
        with open(argv[2], 'r') as in_file:
            lines = in_file.readlines()
            j = 0
            while j < len(lines):
                if lines[j][:4] == 'def ' and lines[j - 1] != '@profile\n':
                    lines.insert(j, '@profile\n')
                    j += 1
                j += 1
            with open(argv[2], 'w') as out_file:
                out_file.writelines(lines)
    
    elif argv[1] == '--postprocess':
        with open(argv[2], 'r') as in_file:
            lines = in_file.readlines()
        lines = [line for line in lines if '@profile' not in line]
        with open(argv[2], 'w') as out_file:
            out_file.writelines(lines)

    elif argv[1] == '--compile':
        functions = list()
        code_lines = list()
        for j in range(2, len(argv)):
            with open(argv[j]) as benchmark_file:
                lines = benchmark_file.read().splitlines()
            functions_file = list() 
            code_lines_file = list() 
            i = 0
            while i < len(lines):
                if lines[i].startswith('Timer unit: '):
                    timer_unit = float(lines[i].split()[2])
                elif lines[i].startswith('Total time: '):
                    try:
                        time = float(lines[i].split()[2]) / timer_unit
                        file = lines[i + 1].split()[1]
                        name = lines[i + 2].split()[1]
                        functions_file.append([name, time, file])
                        i += 2
                    except (IndexError, ValueError), e:
                        pass 
                else:  
                    try:
                        temp_line = lines[i].split()[:5]
                        for j in range(0, 2):
                            temp_line[j] = int(temp_line[j])
                        for j in range(2, 5):
                            temp_line[j] = float(temp_line[j])
                        temp_line.append(name)
                        code_lines_file.append(temp_line)
                    except (IndexError, ValueError), e:
                        pass
                i += 1
            total_time = 0
            for function in functions_file:
                total_time += function[1]
            for function in functions_file:
                function.insert(2, float(function[1] / total_time))
            functions.extend(functions_file)
            code_lines.extend(code_lines_file)
        functions_dict = dict() 
        for function in functions:
            if function[0] not in functions_dict.keys():
                functions_dict[function[0]] = list([list(),list(),function[3]])
            for i in range(2):
                functions_dict[function[0]][i].append(function[i+1])
        functions_summary = list()
        for function in functions_dict.keys():
            functions_summary.append(list([
                function,
                median(functions_dict[function][0]),
                median(functions_dict[function][1]),
                functions_dict[function][2]
                ]))
        with open('codes3d/codes3d.py', 'r') as codes3d_script:
            codes3d_lines = codes3d_script.read().splitlines()
        profile_count = dict() 
        count = 0
        for i in range(len(codes3d_lines)):
            if '@profile' in codes3d_lines[i]: 
                count += 1
            profile_count[i+1] = count
        for i in range(len(codes3d_lines)):
            codes3d_lines[i] = codes3d_lines[i].strip('\t ')
        code_lines_dict = dict()
        for code_line in code_lines:
            if code_line[0] not in code_lines_dict.keys():
                code_lines_dict[code_line[0]] = list([code_line[1],
                                                      list(), 
                                                      list(), 
                                                      list(), 
                                                      code_line[5], 
                                                      codes3d_lines[code_line[0]-1]])
            for i in range(1,4): 
                code_lines_dict[code_line[0]][i].append(code_line[i+1])
        code_lines_summary = list()
        for code_line in code_lines_dict.keys():
            corrected_code_line = str(int(code_line) - profile_count[code_line])
            code_lines_summary.append(list([
                corrected_code_line,
                code_lines_dict[code_line][0],
                median(code_lines_dict[code_line][1]),
                median(code_lines_dict[code_line][2]),
                median(code_lines_dict[code_line][3]),
                code_lines_dict[code_line][4],
                code_lines_dict[code_line][5]]))
        functions_summary.sort(key=lambda function: function[1], reverse=True)
        code_lines_summary.sort(key=lambda code_line: code_line[2], reverse=True)
        timestamp = strftime("%d%m-%H%M%S")
        with open('benchmarks/functions_' + timestamp + '.csv', 'wb') as function_file:
            function_writer = csv_writer(function_file, delimiter=',')
            function_writer.writerow(['Function', 'Execution Time', '% Total Execution Time', 'File'])
            for function in functions_summary:
                function_writer.writerow(function)
        with open('benchmarks/lines_' + timestamp + '.csv', 'wb') as line_file:
            line_writer = csv_writer(line_file, delimiter=',')
            line_writer.writerow(['Line', 'Times Executed', 'Total Execution Time', 'Time per Execution', '% of Function Execution Time', 'Function', 'Code'])
            for code_line in code_lines_summary:
                        line_writer.writerow(code_line)
        print "RESULTS WRITTEN TO:\nbenchmarks/functions_" + timestamp + ".csv\nbenchmarks/lines_" + timestamp + ".csv"
    else:
        quit()

if __name__ == '__main__':
    main()
