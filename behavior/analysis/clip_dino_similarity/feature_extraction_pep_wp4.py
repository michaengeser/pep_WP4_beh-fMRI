import torch
from PIL import Image
from transformers import AutoProcessor, CLIPModel, AutoImageProcessor, AutoModel
import faiss
import os
import numpy as np
import pandas as pd

device = torch.device('cuda' if torch.cuda.is_available() else "cpu")
current_directory = os.getcwd()

# Load CLIP model and processor
processor_clip = AutoProcessor.from_pretrained("openai/clip-vit-base-patch32")
model_clip = CLIPModel.from_pretrained("openai/clip-vit-base-patch32").to(device)

# Load DINOv2 model and processor
processor_dino = AutoImageProcessor.from_pretrained('facebook/dinov2-base')
model_dino = AutoModel.from_pretrained('facebook/dinov2-base').to(device)


# Define a function that normalizes embeddings and add them to the index
def add_vector_to_index(embedding, index):
    # convert embedding to numpy
    vector = embedding.detach().cpu().numpy()
    # Convert to float32 numpy
    vector = np.float32(vector)
    # Normalize vector: important to avoid wrong results when searching
    faiss.normalize_L2(vector)
    # Add to index
    index.add(vector)


def extract_features_clip(image):
    with torch.no_grad():
        inputs = processor_clip(images=image, return_tensors="pt").to(device)
        image_features = model_clip.get_image_features(**inputs)
        return image_features


def extract_features_dino(image):
    with torch.no_grad():
        inputs = processor_dino(images=image, return_tensors="pt").to(device)
        outputs = model_dino(**inputs)
        image_features = outputs.last_hidden_state
        return image_features.mean(dim=1)


# Initialize image types and drawing folders
image_types = ["own", "control"]
drawing_folders = ["draw3D_images_exp", "used_drawings_exp"]
exps = ["1", "2"]

# Base path for the folders
base_path = '../../image_similarities/'

# Loop through image types
for image_type in image_types:

    # Loop through drawing folders
    for drawing_folder in drawing_folders:

        for exp in exps:
            # Initialize tables for extracted features
            feature_table_clip = pd.DataFrame()
            feature_table_dino = pd.DataFrame()

            # Create folder path
            folder_path = os.path.join(current_directory, '..', '..',
            'image_similarities', 'drawings_draw3D', drawing_folder + exp, image_type)

            # Retrieve all filenames
            images = []
            for root, dirs, files in os.walk(folder_path):
                for file in files:
                    if file.endswith('jpg') or file.endswith('png'):
                        images.append(os.path.join(root, file))

            # Print the number of loaded images and the current iteration context
            print(f'Loaded {len(images)} images for {drawing_folder} and {image_type}')

            # for root, dirs, files in os.walk('./' + image_type + '/'):

            # Iterate over the dataset to extract features X2 and store features in tab
            for image_path in images:
                img = Image.open(image_path).convert('L')
                clip_features = extract_features_clip(img)
                dino_features = extract_features_dino(img)
                # store in table
                filename = os.path.basename(image_path)
                filename = os.path.splitext(filename)[0]
                feature_table_clip[filename] = clip_features[0]
                feature_table_dino[filename] = dino_features[0]
                print('Extracted features for', filename)

            # Save the DataFrame to a CSV file
            clip_save_name = image_type + '_' + drawing_folder + exp + "_feature_table_clip.csv"
            feature_table_clip.to_csv(clip_save_name, index=False)
            dino_save_name = image_type + '_' + drawing_folder + exp + "_feature_table_dino.csv"
            feature_table_dino.to_csv(dino_save_name, index=False)

