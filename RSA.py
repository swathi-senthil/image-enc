import numpy as np
from PIL import Image
from Crypto.Cipher import AES, PKCS1_OAEP
from Crypto.PublicKey import RSA
from Crypto.Random import get_random_bytes
from scipy.integrate import odeint

# Define the hyperchaotic Lorenz system
def hyperchaotic_lorenz_system(state, t, sigma, rho, beta, gamma):
    x, y, z = state
    dx_dt = sigma * (y - x)
    dy_dt = x * (rho - z) + y
    dz_dt = x * y - beta * z + gamma * (x ** 2 + y ** 2)
    return [dx_dt, dy_dt, dz_dt]

# Generate chaotic key using the Lorenz system
def generate_hyperchaotic_key(sigma=10, rho=25, beta=8/3, gamma=0.1, timesteps=10000):
    initial_state = [1.0, 1.0, 1.0]
    t = np.arange(0, timesteps * 0.01, 0.01)
    solution = odeint(hyperchaotic_lorenz_system, initial_state, t, args=(sigma, rho, beta, gamma))
    key = np.array(solution[:, 0])  # Using only the x-coordinate as the key
    return key

# RSA key generation
def generate_rsa_key():
    key = RSA.generate(2048)
    public_key = key.publickey().export_key()
    private_key = key.export_key()
    return public_key, private_key

# Circular shift function for preprocessing
def circular_shift(array, shift_amount):
    return np.roll(array, shift_amount)

# Encrypt image using RSA and AES
def encrypt_image(image_path, rsa_public_key, shift_repetitions=100):
    # Load image
    image = Image.open(image_path)
    pixel_array = np.array(image).tobytes()

    # Perform circular shifting as preprocessing
    pixel_array = circular_shift(pixel_array, shift_repetitions)

    # Generate AES key and initialization vector
    aes_key = get_random_bytes(16)
    iv = get_random_bytes(16)

    # Encrypt image pixel values using AES
    cipher_aes = AES.new(aes_key, AES.MODE_CBC, iv)
    encrypted_pixel_values = cipher_aes.encrypt(bytes(pixel_array))

    # Encrypt AES key with RSA
    rsa_public_key_obj = RSA.import_key(rsa_public_key)
    cipher_rsa = PKCS1_OAEP.new(rsa_public_key_obj)
    encrypted_aes_key = cipher_rsa.encrypt(aes_key)

    return encrypted_aes_key, iv, encrypted_pixel_values

# Decrypt image using RSA and AES
def decrypt_image(encrypted_aes_key, iv, encrypted_pixel_values, rsa_private_key, shift_repetitions=100):
    # Decrypt AES key with RSA
    rsa_private_key_obj = RSA.import_key(rsa_private_key)
    cipher_rsa = PKCS1_OAEP.new(rsa_private_key_obj)
    aes_key = cipher_rsa.decrypt(encrypted_aes_key)

    # Decrypt image pixel values using AES
    cipher_aes = AES.new(aes_key, AES.MODE_CBC, iv)
    decrypted_pixel_values = cipher_aes.decrypt(encrypted_pixel_values)

    # Reverse circular shifting for decryption
    decrypted_pixel_values = circular_shift(decrypted_pixel_values, -shift_repetitions)

    return decrypted_pixel_values

# Example usage
if __name__ == "__main__":
    # Generate RSA keys
    rsa_public_key, rsa_private_key = generate_rsa_key()

    # Encrypt image
    encrypted_aes_key, iv, encrypted_pixel_values = encrypt_image('images.jpg', rsa_public_key)

    # Save encrypted image as JPEG
    encrypted_image = Image.frombytes('RGB', (512, 512), encrypted_pixel_values)
    encrypted_image.save('enc2.jpg')

    print("Encrypted image saved as 'enc2.jpg'")

    # Decrypt image
    decrypted_pixel_values = decrypt_image(encrypted_aes_key, iv, encrypted_pixel_values, rsa_private_key)

    # Convert decrypted pixel values to an image
    decrypted_image = Image.frombytes('RGB', (512, 512), decrypted_pixel_values)

    # Save decrypted image as JPEG
    decrypted_image.save('dec2.jpg')
    print("Decrypted image saved as 'dec2.jpg'")