import wandb
import matplotlib.pyplot as plt

api = wandb.Api()
run = api.run("reed-maximilian-university-of-toronto/bioemu-foldedness-finetune/xagoqi7m")

history = run.history()

fig, ax = plt.subplots(figsize=(8, 4))

# Drop NaNs per column before plotting so matplotlib has actual points to connect
train_loss = history[["_epoch", "train/loss"]].dropna()
val_loss = history[["_epoch", "val/loss"]].dropna()

# We'll add a smoothed loss for val loss to make it easier to interpret
less_epochs = val_loss["_epoch"][1:-1]
smoothed_val_loss = (val_loss["val/loss"][2:]+val_loss["val/loss"][1:-1]+val_loss["val/loss"][:-2])/3.

ax.plot(train_loss["_epoch"], train_loss["train/loss"], label="Train loss")
ax.plot(val_loss["_epoch"], val_loss["val/loss"], label="Val loss")
ax.plot(less_epochs, smoothed_val_loss, label="Val loss (Smoothed)", linewidth=4)

ax.set_xlabel("Epoch")
ax.set_ylabel("Loss")
ax.set_title("Foldedness fine-tuning loss")
ax.legend()
ax.axhline(y=0, color="grey", linestyle="--", linewidth=0.8)  # helpful reference for negative loss

fig.savefig("loss_curve.png", dpi=150, bbox_inches="tight")
print("Saved loss_curve.png")
