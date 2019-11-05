from google_drive_downloader import GoogleDriveDownloader as gdd

print('\n\n!!! This may take several minutes !!!\n\n')

gdd.download_file_from_google_drive(file_id='1u9PkZ5UWcKeOaJGjhS-VYm0WIz0vd6cP',
                                    dest_path='./Rlibrary_NeVOmics.zip',
                                    unzip=True)
