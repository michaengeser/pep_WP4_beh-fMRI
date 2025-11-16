from PIL import Image
import os

def resize_images_in_folder(folder_path):
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file.lower().endswith(".jpg"):
                image_path = os.path.join(root, file)
                with Image.open(image_path) as img:
                    # Calculate new size (one quarter)
                    width, height = img.size
                    new_size = (width // 2, height // 2)

                    # Resize the image while maintaining proportions
                    resized_img = img.resize(new_size, Image.LANCZOS)

                    # Save the resized image, overwrite the original
                    resized_img.save(image_path, quality=95)

if __name__ == "__main__":
    # Example usage with relative path from the current directory
    folder_to_process = os.getcwd()

    resize_images_in_folder(folder_to_process)