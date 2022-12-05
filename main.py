import os
from PIL import Image

# List all files in the current directory
files = os.listdir("./pics")

def merge_images(image_one: str, image_two:str, index: int) -> None:
    # Open the images
    image1 = Image.open(image_one)
    image2 = Image.open(image_two)

    # Create a new image with the dimensions of the two images combined
    width = max(image1.width, image2.width)
    height = image1.height + image2.height
    result_image = Image.new('RGB', (width, height))

    # Paste the images into the new image
    result_image.paste(image1, (0, 0))
    result_image.paste(image2, (0, image1.height))

    # Save the new image
    result_image.save('./merged-pics/result-'+str(index)+'.png')

def merge_image(image_one: str, index: int) -> None:
    # Open the images
    image1 = Image.open(image_one)

    # Create a new image with the dimensions of the two images combined
    result_image = Image.new('RGB', (image1.width, image1.height))

    # Paste the images into the new image
    result_image.paste(image1, (0, 0))

    # Save the new image
    result_image.save('./merged-pics/result-'+str(index)+'.png')

# Define the custom sorting function
def sort_strings(string):
  # Extract the numerical part of the string
  number = int(string[:-4])
  
  # Return the numerical part as the sorting key
  return number

# Sort the strings in ascending order using the custom sorting function
sorted_strings = sorted(files, key=sort_strings)
count = 0
for x in sorted_strings:
    if count %2 == 0:
        if count != len(sorted_strings) -1:
            merge_images('./pics/'+sorted_strings[count], './pics/'+sorted_strings[count+1], count)
        else :
            merge_image('./pics/'+sorted_strings[count], count)
    count+=1




print(len(sorted_strings))
