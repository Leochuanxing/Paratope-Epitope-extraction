# -*- coding:utf-8 -*-

from os import listdir
import xlwt


def list_type(path):
    organism_types = {}
    file_list = listdir(path)
    pdb_list = [file for file in file_list if file[-4:] == '.pdb']
    for pdb in pdb_list:
        with open('{}/{}'.format(path, pdb), 'r') as f:
            content = f.readlines()
        stay = False
        organism_scientific, organism_common = '', ''
        for line in content:
            if line[:6] == 'SOURCE':
                if line[11:30] == 'ORGANISM_SCIENTIFIC':
                    organism_scientific = line[31:].strip()
                    if organism_scientific[-1] == ';':
                        organism_scientific = organism_scientific[:-1]
                    else:
                        stay = True
                elif stay:
                    organism_scientific += line[12:].strip()[:-1]
                    stay = False
                if line[11:26] == 'ORGANISM_COMMON':
                    organism_common = line[27:].strip()[:-1]
        if organism_scientific and organism_common:
            organism_type = '{} ({})'.format(organism_scientific, organism_common)
        elif organism_scientific:
            organism_type = organism_scientific
        else:
            organism_type = 'INVALID'
        if organism_type not in organism_types:
            organism_types[organism_type] = 1
            # print(organism_type)
        else:
            organism_types[organism_type] += 1
    return organism_types


def list_type_only_scientific(path):
    organism_types = {}
    file_list = listdir(path)
    pdb_list = [file for file in file_list if file[-4:] == '.pdb']
    for pdb in pdb_list:
        with open('{}/{}'.format(path, pdb), 'r') as f:
            content = f.readlines()
        stay = False
        organism_scientific = '', ''
        for line in content:
            if line[:6] == 'SOURCE':
                if line[11:30] == 'ORGANISM_SCIENTIFIC':
                    organism_scientific = line[31:].strip()
                    if organism_scientific[-1] == ';':
                        organism_scientific = organism_scientific[:-1]
                    else:
                        stay = True
                elif stay:
                    organism_scientific += line[12:].strip()[:-1]
                    stay = False
        if organism_scientific:
            organism_type = organism_scientific
        else:
            organism_type = 'INVALID'
        if organism_type not in organism_types:
            organism_types[organism_type] = 1
            # print(organism_type)
        else:
            organism_types[organism_type] += 1
    return organism_types


def write_result_in_xls(organism_types, out_path='./organism_types.xls'):
    wb = xlwt.Workbook()
    ws = wb.add_sheet('organism_types')
    res_list = []
    for k in organism_types:
        res_list.append((k, organism_types[k]))
    res_list.sort(key=lambda x: x[1], reverse=True)
    for row in range(len(res_list)):
        ws.write(row, 0, res_list[row][0])
        ws.write(row, 1, res_list[row][1])
    wb.save(out_path)


if __name__ == '__main__':
    res = list_type(
        './data/resolusion 4.0 complex True polypeptide > 5/structure')
    print(res)
    write_result_in_xls(res)

    res = list_type_only_scientific(
        './data/resolusion 4.0 complex True polypeptide > 5/structure')
    print(res)
    write_result_in_xls(res, './organism_types_only_scientific.xls')
