import itertools
import csv
import pandas as pd
import random as rnd


def generate_condition_template():
    subjects_numbers = [num for num in range(101, 141+1) if num not in [107, 108, 110, 115, 128, 133]]
    subjects_numbers = [str(num) for num in subjects_numbers]

    # Ensure that there are exactly 36 subjects
    if len(subjects_numbers) != 35:
        raise ValueError("Number of subjects is not equal to 35.")

    # make pairs and add id
    pairs = list(itertools.combinations(subjects_numbers, 2))
    pairs = [(pair[0], pair[1], pair_id) for pair_id, pair in enumerate(pairs, start=1)]

    # Create a CSV file and write header
    with open('condition_template.csv', 'w', newline='') as csvfile:
        fieldnames = ['sub1', 'sub2', 'pairID', 'trialID']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        # Write data to the CSV file
        for pair in pairs:
            writer.writerow({'sub1': pair[0], 'sub2': pair[1], 'pairID': pair[2], 'trialID': pair[2]})

        # Set random seed for reproducibility
        rnd.seed(1)

        trial_id = pairs[-1][2]
        # Add repeats of randomly selected rows
        for _ in range(5):
            random_row = rnd.choice(pairs)
            for _ in range(4):
                trial_id += 1
                writer.writerow({'sub1': random_row[0], 'sub2': random_row[1], 'pairID': random_row[2], 'trialID': trial_id})


def create_condition_files(template_csv, category_name, category_abbreviation):
    # Read the template CSV file
    template_df = pd.read_csv(template_csv)

    # Add 'stim1' and 'stim2' rows
    template_df['stim1'] = ''
    template_df['stim2'] = ''

    # Construct image paths for 'stim1' and 'stim2'
    for index, row in template_df.iterrows():
        stim1_path = f'stimuli/own/{category_name}/{row["sub1"]}_{category_abbreviation}.jpg'
        stim2_path = f'stimuli/own/{category_name}/{row["sub2"]}_{category_abbreviation}.jpg'
        template_df.at[index, 'stim1'] = stim1_path
        template_df.at[index, 'stim2'] = stim2_path

    # Save the files with specified names
    file_name = f'conditions_{category_name}.csv'
    template_df.to_csv(file_name, index=False)


if __name__ == "__main__":
    generate_condition_template()

    template_csv_path = 'condition_template.csv'
    category_name_kitchen = 'kitchen'
    category_name_bathroom = 'bathroom'
    category_abbreviation_kit = 'kit'
    category_abbreviation_bat = 'bat'

    create_condition_files(template_csv_path, category_name_kitchen, category_abbreviation_kit)
    create_condition_files(template_csv_path, category_name_bathroom, category_abbreviation_bat)
