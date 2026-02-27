import wandb
import matplotlib.pyplot as plt

api = wandb.Api()
run = api.run("reed-maximilian/bioemu-foldedness-finetune/pious-spaceship-5")

# Get the full history as a dataframe
history = run.history()

# Plot and save however you like
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(history["train/loss"], label="Train")
ax.plot(history["val/loss"], label="Val")
ax.set_xlabel("Step")
ax.set_ylabel("Loss")
ax.legend()
fig.savefig("loss_curve.png", dpi=150, bbox_inches="tight")
