import wandb
import matplotlib.pyplot as plt

api = wandb.Api()
run = api.run("reed-maximilian-university-of-toronto/bioemu-foldedness-finetune/xagoqi7m")

history = run.history()

fig, ax = plt.subplots(figsize=(8, 4))

# Drop NaNs per column before plotting so matplotlib has actual points to connect
train_loss = history[["_step", "train/loss"]].dropna()
val_loss = history[["_step", "val/loss"]].dropna()

ax.plot(train_loss["_step"], train_loss["train/loss"], label="Train loss")
ax.plot(val_loss["_step"], val_loss["val/loss"], label="Val loss")

ax.set_xlabel("Step")
ax.set_ylabel("Loss")
ax.set_title("Foldedness fine-tuning loss")
ax.legend()
ax.axhline(y=0, color="grey", linestyle="--", linewidth=0.8)  # helpful reference for negative loss

fig.savefig("loss_curve.png", dpi=150, bbox_inches="tight")
print("Saved loss_curve.png")
