import os

cmd = "pipreqs ./pyCPTM/ --force"
os.system(cmd)

# Alternative: pip freeze > requirements.txt