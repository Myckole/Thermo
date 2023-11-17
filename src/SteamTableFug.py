def read_file(file_name):
    with open(file_name, 'r') as file:
        data = file.readlines()
    data = [list(map(float, line.split())) for line in data]
    return data
