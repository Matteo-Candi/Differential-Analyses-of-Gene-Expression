import json

def delete_keys_with_words(data, words_to_delete):
    if isinstance(data, dict):
        for key in list(data.keys()):
            if any(word in key for word in words_to_delete):
                del data[key]
            else:
                delete_keys_with_words(data[key], words_to_delete)
    elif isinstance(data, list):
        for item in data:
            delete_keys_with_words(item, words_to_delete)

# Load the JSON data
with open('biospecimen.project-tcga-chol.2024-04-18.json', 'r') as file:
    json_data = json.load(file)

# Define the list of words to delete
words_to_delete = ['created_datetime',
                    'submitter_id',
                    'updated_datetime',
                    'aliquot_id',
                    'center_id',
                    'name',
                    'portion_id',
                    'pathology_report_uuid',
                    'sample_id',
                    'analyte_id',
                    'slide_id',
                    'case_id',
                    'annotation_id',
                    'entity_id'
                    ]

# Delete keys containing specified words
delete_keys_with_words(json_data, words_to_delete)

# Write the modified data back to the file
with open('filtered_biospecimen.json', 'w') as file:
    json.dump(json_data, file, indent=4)
