import csv

class Parameters:
    def __init__(self):
        self.parameters = {}

    # Function to add parameters as object instances 
    def add_parameter(self, name, value, min_bound = None, max_bound = None):
        value = self._apply_bounds(value, min_bound, max_bound)
        self.parameters[name] = {
            'value': value,
            'min_bound': min_bound,
            'max_bound': max_bound
        }

    # Function to apply the bounds on the parameters 
    def _apply_bounds(self, value, min_bound, max_bound):
        if min_bound is not None and value < min_bound:
            return min_bound
        if max_bound is not None and value > max_bound:
            return max_bound
        return value

    # Function to export parameters as a dictionary 
    def to_dict(self):
        return {name: param['value'] for name, param in self.parameters.items()}

    # Function to export parameters as a csv file 
    def export_to_csv(self, filename):
        with open(filename, mode = 'w', newline = '') as file:
            writer = csv.writer(file)
            writer.writerow(['name', 'value', 'min_bound', 'max_bound'])
            for name, param in self.parameters.items():
                writer.writerow([name, param['value'], param['min_bound'], param['max_bound']])

    # Function to import parameters from a csv file 
    @staticmethod
    def import_from_csv(filename):
        params = Parameters()
        with open(filename, mode = 'r') as file:
            reader = csv.DictReader(file)
            for row in reader:
                name = row['name']
                value = float(row['value'])  # Assuming values are numeric
                min_bound = float(row['min_bound']) if row['min_bound'] else None
                max_bound = float(row['max_bound']) if row['max_bound'] else None
                params.add_parameter(name, value, min_bound, max_bound)
        return params